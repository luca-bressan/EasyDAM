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
            profile = zeros(24)
            profile[bid.market_period]=bid.quantity*(bid.purpose==EasyDAM.sell ? 1 : -1)
            push!(flattened_bids,(id=bid.id, price=bid.price, MAR=0.0, zone=bid.zone, profile=profile))
        end
        if isa(bid,BlockBid)
            profile = zeros(24)
            for slice in bid.block
                profile[slice.market_period]=slice.quantity*(bid.purpose==EasyDAM.sell ? 1 : -1)
            end
            push!(flattened_bids,(id=bid.id, price=bid.price, MAR=bid.MAR, zone=bid.zone, profile=profile))
        end
    end
    return flattened_bids
end

# Market building methods
function get_zonal_limits(lines::Vector{Line}, zones::Vector{String})::DataFrames.DataFrame
    export_profiles = Vector{Vector{Float64}}()
    import_profiles = Vector{Vector{Float64}}()

    # Process each zone
    for zone in zones
        export_profile = zeros(24)  # Assuming 24 hours as per the original logic
        import_profile = zeros(24)
        
        # Aggregate profiles for each zone
        for line in lines
            if line.source == zone
                export_profile[line.market_period] += line.ATC
            elseif line.destination == zone
                import_profile[line.market_period] -= line.ATC
            end
        end

        push!(export_profiles, export_profile)
        push!(import_profiles, import_profile)
    end

    # Create DataFrame
    return DataFrames.DataFrame(
        zone = zones,
        export_limit = export_profiles,
        import_limit = import_profiles
    )
end

function get_line_profiles(lines::Vector{Line}, zones::Vector{String})::DataFrames.DataFrame
    upstream_profiles = Vector{Vector{Float64}}()
    downstream_profiles = Vector{Vector{Float64}}()
    line_ids = Set{Tuple{String, String}}()
    line_froms = Vector{String}()
    line_tos = Vector{String}()
    ids = Vector{String}()

    # Populate unique line_ids
    for line in lines
        line_tuple = (line.source, line.destination)
        reverse_tuple = (line.destination, line.source)
        if !(line_tuple in line_ids || reverse_tuple in line_ids)
            push!(line_ids, line_tuple)
        end
    end

    # Process each unique line_id
    for (line_from, line_to) in line_ids
        upstream_profile = zeros(24)  # Assuming 24 hours as per the original logic
        downstream_profile = zeros(24)
        
        # Aggregate profiles
        for line in lines
            if line.source == line_from && line.destination == line_to
                upstream_profile[line.market_period] += line.ATC
            elseif line.source == line_to && line.destination == line_from
                downstream_profile[line.market_period] += line.ATC
            end
        end

        push!(line_froms, line_from)
        push!(line_tos, line_to)
        push!(ids, line_from * line_to)
        push!(upstream_profiles, upstream_profile)
        push!(downstream_profiles, downstream_profile)
    end

    # Create DataFrame
    return DataFrames.DataFrame(
        line_id = ids,
        line_from = line_froms,
        line_to = line_tos,
        upstream_profile = upstream_profiles,
        downstream_profile = downstream_profiles
    )
end
# Home of the magical gnomes
function build_and_solve_market(bids::DataFrames.DataFrame,zonal_limits::DataFrames.DataFrame,line_profiles::DataFrames.DataFrame)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(model, x[bids.id] >= 0.0)
    for bid in eachrow(bids)
        JuMP.@constraint(model, x[bid.id] <= 1.0) #TODO: add MAR usage
    end
    JuMP.@variable(model, s[bids.id] >= 0.0)
    bids.x = Array(x)
    bids.s = Array(s)
    JuMP.@variable(model, p[i=1:24,zonal_limits.zone])
    zonal_limits.p = Array([p[1:24, zl.zone] for zl in eachrow(zonal_limits)])
    for zl in eachrow(zonal_limits)
       JuMP.@constraint(model, p[:,zl.zone] .<= 9999.0)
       JuMP.@constraint(model, p[:,zl.zone] .>= -9999.0)
    end
    market = DataFrames.leftjoin(bids,zonal_limits,on=:zone)
    JuMP.@variable(model, upstream_flow[i=1:24,line_profiles.line_id])
    JuMP.@variable(model, downstream_flow[i=1:24,line_profiles.line_id])
    JuMP.@variable(model, upstream_profit[i=1:24,line_profiles.line_id])
    JuMP.@variable(model, downstream_profit[i=1:24,line_profiles.line_id])
    for line in eachrow(line_profiles)
        JuMP.@constraint(model, upstream_profit[:,line.line_id] .>= zeros(24))
        JuMP.@constraint(model, downstream_profit[:,line.line_id] .>= zeros(24))
        JuMP.@constraint(model, upstream_flow[:,line.line_id] .>= zeros(24))
        JuMP.@constraint(model, downstream_flow[:,line.line_id] .>= zeros(24))
        JuMP.@constraint(model, upstream_flow[:,line.line_id] .<= line.upstream_profile)
        JuMP.@constraint(model, downstream_flow[:,line.line_id] .<= line.downstream_profile)
    end
    JuMP.@expression(model, flow[i=1:24,line=line_profiles.line_id], upstream_flow[i,line] .- downstream_flow[i,line])
    
    
    # grid_balance_constraint = JuMP.@constraint(model,
    #     grid_balance_constraint,
    #     sum(market.x .* market.profile, dims=1)[1] .== zeros(24))
    
    xprt = Dict{String, Vector{JuMP.AffExpr}}()
    for zl in eachrow(zonal_limits)
        zone = zl.zone  
        xprt[zone] = Vector{JuMP.AffExpr}(undef, 24)  
        for i in eachindex(xprt[zone])
            xprt[zone][i] = JuMP.AffExpr(0.0)
        end
        for line in eachrow(line_profiles)
            if line.line_from == zone
                xprt[zone] = xprt[zone] .+ flow[:, line.line_id]
            elseif line.line_to == zone
                xprt[zone] = xprt[zone] .- flow[:, line.line_id]
            end
        end
    end

    for zl in eachrow(zonal_limits)
        # println(zl.zone * "-" * string(isempty(market[market.zone.==zl.zone,:])))
        if !isempty(market[market.zone.==zl.zone,:])
            JuMP.@constraint(model,
            sum(market[market.zone.==zl.zone,:].x .* market[market.zone.==zl.zone,:].profile, dims=1)[1] .>= xprt[zl.zone])
            JuMP.@constraint(model,
            sum(market[market.zone.==zl.zone,:].x .* market[market.zone.==zl.zone,:].profile, dims=1)[1] .<= xprt[zl.zone])
        else
            JuMP.@constraint(model,
            xprt[zl.zone] .== zeros(24))
        end
    end
    JuMP.@constraint(model,
        [b in eachrow(market)],
        b.s >= sum(b.profile .* (b.p .- b.price))
    )

    JuMP.@constraint(model,
        [l in eachrow(line_profiles)],
        p[:,l.line_to] .- p[:,l.line_from] .<= upstream_profit[:,l.line_id]
    )

    JuMP.@constraint(model,
        [l in eachrow(line_profiles)],
        p[:,l.line_from] .- p[:,l.line_to] .<= downstream_profit[:,l.line_id]
    )

    JuMP.@constraint(
        model,
        duality_constraint,
        sum(market.s)
        + sum([sum(upstream_profit[:, line.line_id] .* line.upstream_profile) for line in eachrow(line_profiles)])
        + sum([sum(downstream_profit[:, line.line_id] .* line.downstream_profile) for line in eachrow(line_profiles)])
        == -sum(sum(market.x .* market.profile .* market.price))
    )

    JuMP.@objective(model, Max, -sum(sum(market.x .* market.profile .* market.price)))
    JuMP.optimize!(model)
    @assert JuMP.is_solved_and_feasible(model)
    println(JuMP.solution_summary(model))
    for l in eachrow(line_profiles)
        if l.line_id == "MALTSICI"
            println(l.line_id * "- flow: " *string(JuMP.value.(flow[:,l.line_id]))* "- uflow: " *string(JuMP.value.(upstream_flow[:,l.line_id]))* "- dflow: " *string(JuMP.value.(downstream_flow[:,l.line_id])))
            println(l.line_id * "- uflow: " *string(JuMP.value.(upstream_profit[:,l.line_id]))* "- dflow: " *string(JuMP.value.(downstream_profit[:,l.line_id])))
        end
    end
    market.zp = [Array(JuMP.value.(m.p)) for m in eachrow(market)]
    #market.mcp_unc = [Array(JuMP.value.(MCP_UNC[1:24])) for m in eachrow(market)]
    market.x = [JuMP.value(m.x) for m in eachrow(market)]
    market.s = [JuMP.value(m.s) for m in eachrow(market)]
    zonal_limits.zp = [Array(JuMP.value.(zl.p)) for zl in eachrow(zonal_limits)]
    return (zonal_limits[:,[:zone,:zp]],market[:, [:zone,:profile,:price,:x,:s,:zp]],model)
    #return (market[:, [:zone,:profile,:x,:p]],model)
    
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
