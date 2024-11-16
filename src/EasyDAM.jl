module EasyDAM
import YAML
import StructTypes
import DataFrames.DataFrames
import JuMP
import HiGHS

# Type declarations
struct Line
    source::String
    destination::String
    ATC::Float64
    market_period::UInt8
end
StructTypes.StructType(::Type{Line}) = StructTypes.Struct()

@enum Purpose buy sell
StructTypes.StructType(::Type{Purpose}) = StructTypes.StringType()

@enum BidType block simple
StructTypes.StructType(::Type{BidType}) = StructTypes.StringType()

struct BlockSlice
    market_period::UInt64
    quantity::Float64
end
StructTypes.StructType(::Type{BlockSlice}) = StructTypes.Struct()

abstract type Bid end

struct BlockBid <: Bid
    id::UInt64
    block::Vector{BlockSlice}
    price::Float64
    type::EasyDAM.BidType
    purpose::EasyDAM.Purpose
    zone::String
    MAR::Float64
end
StructTypes.StructType(::Type{BlockBid}) = StructTypes.Struct()

struct SimpleBid <: Bid
    id::UInt64
    quantity::Float64
    price::Float64
    market_period::UInt8
    type::EasyDAM.BidType
    purpose::EasyDAM.Purpose
    zone::String
end
StructTypes.StructType(::Type{SimpleBid}) = StructTypes.Struct()

# Exception declarations
struct DuplicateLineError <: Exception
    msg::String
end
struct DuplicateBidError <: Exception
    msg::String
end
struct BadFileError <: Exception
    msg::String
end

# File parsing methods
function load_lines(filepath::String)::Vector{Line}
    if !isfile(filepath) || !(splitext(filepath)[end] in [".yaml",".yml"])
        throw(BadFileError(filepath * ": not a valid file."))
    end
    lines_dict = YAML.load_file(filepath; dicttype=Dict{Symbol,Any})
    lines = Line[]
    if !(:lines in keys(lines_dict))
        throw(BadFileError(filepath * ": this file might be ill-formed."))
    end
    for line in lines_dict[:lines]
        try
            line = StructTypes.constructfrom(Line,line)
            if (line.source,line.destination,line.market_period) in [(l.source,l.destination,l.market_period) for l in lines]
                throw(DuplicateLineError("A line originating in " * line.source *
                " and terminating in " *
                line.destination *
                " for market period " *
                string(line.market_period) *
                " has been specified twice in " *
                filepath * "."))
            end
            push!(lines,line)
        catch e
            if isa(e,MethodError)
                throw(BadFileError(filepath * ": this file might be ill-formed."))
            else
                rethrow(e)
            end
        end
    end
    return lines
end

function load_bids(filepath::String)::Vector{Bid}
    if !isfile(filepath) || !(splitext(filepath)[end] in [".yaml", ".yml"])
        throw(BadFileError(filepath * ": not a valid file."))
    end
    
    bids_dict = YAML.load_file(filepath; dicttype=Dict{Symbol, Any})
    simple_bids = SimpleBid[]
    block_bids = BlockBid[]
    ids = UInt64[]
    if !(:bids in keys(bids_dict))
        throw(BadFileError(filepath * ": this file might be ill-formed."))
    end
    for bid in bids_dict[:bids]
        try  # attempt to parse as a block bid...
            block_bid = StructTypes.constructfrom(BlockBid, bid)
            if block_bid.id in ids
                throw(DuplicateBidError("A bid with id number " *
                                        string(block_bid.id) *
                                        " has been specified twice in " *
                                        filepath * "."))
            end
            push!(block_bids, block_bid)
            push!(ids, block_bid.id)
        catch e
            if isa(e, MethodError)  # ... it failed! Attempt to parse as a simple bid
                try
                    simple_bid = StructTypes.constructfrom(SimpleBid, bid)
                    if simple_bid.id in ids
                        throw(DuplicateBidError("A bid with id number " *
                                                string(simple_bid.id) *
                                                " has been specified twice in " *
                                                filepath * "."))
                    end
                    push!(simple_bids, simple_bid)
                    push!(ids, simple_bid.id)
                catch f
                    if isa(f, MethodError)
                        throw(BadFileError(filepath * ": this file might be ill-formed."))
                    else
                        rethrow(f)
                    end
                end 
            else
                rethrow(e)
            end
        end 
    end

    return [simple_bids; block_bids]
end

function flatten_bids(bids::Vector{Bid})::DataFrames.DataFrame
    flattened_bids = DataFrames.DataFrame(id=UInt64[],price=Float64[],MAR=Float64[],zone=String[],profile=Vector{Float64}[])
    for bid in bids
        if isa(bid,SimpleBid)
            profile = zeros(100)
            profile[bid.market_period]=bid.quantity*(bid.purpose==EasyDAM.sell ? 1 : -1)
            push!(flattened_bids,(id=bid.id, price=bid.price, MAR=0.0, zone=bid.zone, profile=profile))
        end
        if isa(bid,BlockBid)
            profile = zeros(100)
            for slice in bid.block
                profile[slice.market_period]=slice.quantity*(bid.purpose==EasyDAM.sell ? 1 : -1)
            end
            push!(flattened_bids,(id=bid.id, price=bid.price, MAR=bid.MAR, zone=bid.zone, profile=profile))
        end
    end
    return flattened_bids
end

# Market building methods
function get_zonal_limits(lines::Vector{Line},zones::Vector{String})::DataFrames.DataFrame
    export_profiles = Vector{Float64}[]
    import_profiles = Vector{Float64}[]
    for zone in zones
        export_profile = zeros(100)
        import_profile = zeros(100)
        for line in lines
            if line.source == zone
                export_profile[line.market_period] += line.ATC
            end
            if line.destination == zone
                import_profile[line.market_period] -= line.ATC
            end
        end
        push!(export_profiles,export_profile)
        push!(import_profiles,import_profile)
    end
    return DataFrames.DataFrame(zone=zones,export_limit=export_profiles,import_limit=import_profiles)
end

# Home of the magical gnomes
function build_and_solve_market(bids::DataFrames.DataFrame,zonal_limits::DataFrames.DataFrame)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(model, x[bids.id])
    for bid in eachrow(bids)
        JuMP.@constraint(model, x[bid.id] in JuMP.Semicontinuous(bid.MAR, 1.0))
    end
    JuMP.@variable(model, s[bids.id] >= 0)
    bids.x = Array(x)
    bids.s = Array(s)
    JuMP.@variable(model, MCP_UNC[i=1:100])
    JuMP.@variable(model, u[i=1:100,zonal_limits.zone] >=0)
    JuMP.@variable(model, v[i=1:100,zonal_limits.zone] >=0)
    zonal_limits.p = Array([v[1:100, zl.zone].-u[1:100, zl.zone].+MCP_UNC[1:100] for zl in eachrow(zonal_limits)])
    zonal_limits.cc = Array([v[1:100, zl.zone].-u[1:100, zl.zone] for zl in eachrow(zonal_limits)])
    zonal_limits.u = Array([u[1:100, zl.zone] for zl in eachrow(zonal_limits)])
    zonal_limits.v = Array([v[1:100, zl.zone] for zl in eachrow(zonal_limits)])
    zonal_limits.MCP_UNC = Array([MCP_UNC[1:100] for zl in eachrow(zonal_limits)])
    zonal_limits2 = zonal_limits
    filter!(zl -> zl.zone in bids.zone, zonal_limits2)
    market = DataFrames.leftjoin(bids,zonal_limits,on=:zone)
    sourcing_constraint = JuMP.@constraint(model,
        [zl in eachrow(zonal_limits2)],
        sum(market[market.zone.==zl.zone,:].x .* market[market.zone.==zl.zone,:].profile, dims=1)[1] .<= zl.export_limit)
    sinking_constraint = JuMP.@constraint(model,
        [zl in eachrow(zonal_limits2)],
        sum(market[market.zone.==zl.zone,:].x .* market[market.zone.==zl.zone,:].profile, dims=1)[1] .>= zl.import_limit)
    grid_balance_constraint = JuMP.@constraint(model,
        grid_balance_constraint,
        sum(market.x .* market.profile, dims=1)[1] .== zeros(100))
    surplus_constraint = JuMP.@constraint(model,
        [b in eachrow(market)],
        b.s >= sum(b.profile .* (b.p .- b.price))
    )
    JuMP.@constraint(model,
        duality_constraint,
        sum(market.s) == sum([sum(m.x .* m.profile .* m.price) for m in eachrow(market)]) 
                        - sum([sum(zl.import_limit .* zl.u) for zl in eachrow(zonal_limits)])
                        + sum([sum(zl.export_limit .* zl.v) for zl in eachrow(zonal_limits)])
    )

    JuMP.@objective(model, Max, -sum(sum(market.x .* market.profile .* market.price)))
    JuMP.optimize!(model)
    @assert JuMP.is_solved_and_feasible(model)
    JuMP.solution_summary(model)
    market.p = [Array(JuMP.value.(m.p)) for m in eachrow(market)]
    market.x = [JuMP.value(m.x) for m in eachrow(market)]
    market.s = [JuMP.value(m.s) for m in eachrow(market)]
    zonal_limits.p = [Array(JuMP.value.(zl.p)) for zl in eachrow(zonal_limits)]
    zonal_limits.cc = [Array(JuMP.value.(zl.cc)) for zl in eachrow(zonal_limits)]
    zonal_limits.MCP_UNC = [Array(JuMP.value.(zl.MCP_UNC)) for zl in eachrow(zonal_limits)]
    return (zonal_limits[:, [:zone,:MCP_UNC,:cc,:p]],market[:, [:zone,:profile,:price,:x,:p,:s]])
    
end

function recover_transfers(market_results::DataFrames.DataFrame)::DataFrames.DataFrame
end

# File output
function publish_results(market_results::DataFrames.DataFrame, transfers::DataFrames.DataFrame, filepath::String)
end

export load_lines
export load_bids
export Line
export BadFileError
export DuplicateLineError

end # module EasyDAM
