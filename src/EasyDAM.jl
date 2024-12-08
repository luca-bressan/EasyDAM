module EasyDAM
using YAML
using StructTypes
using DataFrames
using JuMP
using HiGHS

# Global config
global const NUM_TIMESTEPS::UInt64 = 24

# Type declarations
struct Line
    source::String
    destination::String
    ATC::Float64
    market_period::UInt64
end
StructType(::Type{Line}) = Struct()

@enum Purpose buy sell
StructType(::Type{Purpose}) = StringType()

@enum BidType block simple
StructType(::Type{BidType}) = StringType()

struct BlockSlice
    market_period::UInt64
    quantity::Float64
end
StructType(::Type{BlockSlice}) = Struct()

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
StructType(::Type{BlockBid}) = Struct()

struct SimpleBid <: Bid
    id::UInt64
    quantity::Float64
    price::Float64
    market_period::UInt64
    type::EasyDAM.BidType
    purpose::EasyDAM.Purpose
    zone::String
end
StructType(::Type{SimpleBid}) = Struct()

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

"""
    load_lines(filepath::String) -> Vector{Line}

Loads and validates lines from a YAML file. Ensures no duplicate lines and 
that the file is well-formed.

# Arguments
- `filepath::String`: Path to the YAML file containing the lines.

# Returns
- A vector of `Line` objects parsed from the file.

# Throws
- `BadFileError`: If the file is invalid or ill-formed.
- `DuplicateLineError`: If duplicate lines are found in the file.
"""
function load_lines(filepath::String)::Vector{Line}
    if !isfile(filepath) || !(splitext(filepath)[end] in [".yaml", ".yml"])
        throw(BadFileError("$filepath: not a valid file."))
    end

    lines_dict = YAML.load_file(filepath; dicttype=Dict{Symbol, Any})
    if !haskey(lines_dict, :lines)
        throw(BadFileError("$filepath: this file might be ill-formed."))
    end

    lines = Line[]
    for raw_line in lines_dict[:lines]
        line = try
            StructTypes.constructfrom(Line, raw_line)
        catch e
            if isa(e, MethodError)
                throw(BadFileError("$filepath: this file might be ill-formed."))
            else
                rethrow(e)
            end
        end

        if is_duplicate_line(line, lines)
            throw(DuplicateLineError(
                "A line originating in $(line.source) and terminating in \
                $(line.destination) for market period $(line.market_period) \
                has been specified twice in $filepath."
            ))
        end

        push!(lines, line)
    end

    return lines
end

function is_duplicate_line(new_line::Line, existing_lines::Vector{Line})::Bool
    return (new_line.source, new_line.destination, new_line.market_period) in 
           [(line.source, line.destination, line.market_period) for line in 
           existing_lines]
end

"""
    load_bids(filepath::String) -> Vector{Bid}

Loads and validates bids from a YAML file. Ensures no duplicate bids and that \
the file is well-formed.

# Arguments
- `filepath::String`: Path to the YAML file containing the bids.

# Returns
- A vector of `Bid` objects (both `SimpleBid` and `BlockBid`) parsed from the \
  file.

# Throws
- `BadFileError`: If the file is invalid or ill-formed.
- `DuplicateBidError`: If duplicate bids are found in the file.
"""
function load_bids(filepath::String)::Vector{Bid}
    if !isfile(filepath) || !(splitext(filepath)[end] in [".yaml", ".yml"])
        throw(BadFileError("$filepath: not a valid file."))
    end

    bids_dict = YAML.load_file(filepath; dicttype=Dict{Symbol, Any})
    simple_bids = SimpleBid[]
    block_bids = BlockBid[]
    ids = UInt64[]
    if !haskey(bids_dict, :bids)
        throw(BadFileError("$filepath: this file might be ill-formed."))
    end

    for bid in bids_dict[:bids]
        try
            block_bid = StructTypes.constructfrom(BlockBid, bid)
            if block_bid.id in ids
                throw(DuplicateBidError(
                    "A bid with id number $(block_bid.id) has been specified \
                    twice in $filepath."
                ))
            end
            push!(block_bids, block_bid)
            push!(ids, block_bid.id)
        catch e
            if isa(e, MethodError)
                try
                    simple_bid = StructTypes.constructfrom(SimpleBid, bid)
                    if simple_bid.id in ids
                        throw(DuplicateBidError(
                            "A bid with id number $(simple_bid.id) has been \
                            specified twice in $filepath."
                        ))
                    end
                    push!(simple_bids, simple_bid)
                    push!(ids, simple_bid.id)
                catch f
                    if isa(f, MethodError)
                        throw(BadFileError(
                            "$filepath: this file might be ill-formed."
                        ))
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

"""
    flatten_bids(bids::Vector{Bid}) -> DataFrame

Converts a vector of `Bid` objects into a flattened `DataFrame` format with \
fields such as `id`, `price`, `MAR`, `zone`, and `profile`.

# Arguments
- `bids::Vector{Bid}`: The vector of `Bid` objects to flatten.

# Returns
- A `DataFrame` with the flattened representation of the bids.
"""
function flatten_bids(bids::Vector{Bid})::DataFrame
    flattened_bids = DataFrame(
        id=UInt64[], price=Float64[], MAR=Float64[], zone=String[], 
        profile=Vector{Float64}[]
    )

    for bid in bids
        if isa(bid, SimpleBid)
            profile = zeros(NUM_TIMESTEPS)
            profile[bid.market_period] = bid.quantity *
                (bid.purpose == EasyDAM.sell ? 1 : -1)
            push!(flattened_bids, (
                id=bid.id, price=bid.price, MAR=0.0, zone=bid.zone, 
                profile=profile
            ))
        end

        if isa(bid, BlockBid)
            profile = zeros(NUM_TIMESTEPS)
            for slice in bid.block
                profile[slice.market_period] = slice.quantity *
                    (bid.purpose == EasyDAM.sell ? 1 : -1)
            end
            push!(flattened_bids, (
                id=bid.id, price=bid.price, MAR=bid.MAR, zone=bid.zone, 
                profile=profile
            ))
        end
    end

    return flattened_bids
end

# Market building methods
"""
    get_zonal_limits(lines::Vector{Line}, zones::Vector{String}) -> DataFrame

Calculates the zonal export and import limits based on the given lines and
zones.

# Arguments
- `lines::Vector{Line}`: A vector of `Line` objects defining the network.
- `zones::Vector{String}`: A vector of zone names to calculate limits for.

# Returns
- A `DataFrame` containing zones, export limits, and import limits.
"""
function get_zonal_limits(lines::Vector{Line}, zones::Vector{String})::DataFrame
    export_profiles = Vector{Vector{Float64}}()
    import_profiles = Vector{Vector{Float64}}()

    for zone in zones
        export_profile = zeros(NUM_TIMESTEPS)
        import_profile = zeros(NUM_TIMESTEPS)

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

    return DataFrame(
        zone = zones,
        export_limit = export_profiles,
        import_limit = import_profiles
    )
end

"""
    get_line_profiles(lines::Vector{Line}, zones::Vector{String}) -> DataFrame

Calculates upstream and downstream profiles for each unique line defined by
its source and destination.

# Arguments
- `lines::Vector{Line}`: A vector of `Line` objects defining the network.
- `zones::Vector{String}`: A vector of zone names (though zones are not used in
  this function).

# Returns
- A `DataFrame` containing line IDs, sources, destinations, upstream profiles,
  and downstream profiles.
"""
function get_line_profiles(lines::Vector{Line}, zones::Vector{String})::DataFrame
    upstream_profiles = Vector{Vector{Float64}}()
    downstream_profiles = Vector{Vector{Float64}}()
    line_ids = Set{Tuple{String, String}}()
    line_froms = Vector{String}()
    line_tos = Vector{String}()
    ids = Vector{String}()

    for line in lines
        line_tuple = (line.source, line.destination)
        reverse_tuple = (line.destination, line.source)
        if !(line_tuple in line_ids || reverse_tuple in line_ids)
            push!(line_ids, line_tuple)
        end
    end

    for (line_from, line_to) in line_ids
        upstream_profile = zeros(NUM_TIMESTEPS)
        downstream_profile = zeros(NUM_TIMESTEPS)

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

    return DataFrame(
        line_id = ids,
        line_from = line_froms,
        line_to = line_tos,
        upstream_profile = upstream_profiles,
        downstream_profile = downstream_profiles
    )
end

# Home of the magical gnomes
function build_and_solve_market(bids::DataFrame,zonal_limits::DataFrame,line_profiles::DataFrame)::Tuple{DataFrame,DataFrame,DataFrame}
    model = Model(HiGHS.Optimizer)
    @variable(model, x[bids.id] >= 0.0)
    for bid in eachrow(bids)
        @constraint(model, x[bid.id] <= 1.0) #TODO: add MAR usage
    end
    @variable(model, s[bids.id] >= 0.0)
    bids.x = Array(x)
    bids.s = Array(s)
    @variable(model, p[i=1:NUM_TIMESTEPS,zonal_limits.zone])
    zonal_limits.p = Array([p[:, zl.zone] for zl in eachrow(zonal_limits)])
    for zl in eachrow(zonal_limits)
       @constraint(model, p[:,zl.zone] .<= 9999.0)
       @constraint(model, p[:,zl.zone] .>= -9999.0)
    end
    market = leftjoin(bids,zonal_limits,on=:zone)
    @variable(model, upstream_flow[i=1:NUM_TIMESTEPS,line_profiles.line_id])
    @variable(model, downstream_flow[i=1:NUM_TIMESTEPS,line_profiles.line_id])
    @variable(model, upstream_profit[i=1:NUM_TIMESTEPS,line_profiles.line_id])
    @variable(model, downstream_profit[i=1:NUM_TIMESTEPS,line_profiles.line_id])
    for line in eachrow(line_profiles)
        @constraint(model, upstream_profit[:,line.line_id] .>= zeros(NUM_TIMESTEPS))
        @constraint(model, downstream_profit[:,line.line_id] .>= zeros(NUM_TIMESTEPS))
        @constraint(model, upstream_flow[:,line.line_id] .>= zeros(NUM_TIMESTEPS))
        @constraint(model, downstream_flow[:,line.line_id] .>= zeros(NUM_TIMESTEPS))
        @constraint(model, upstream_flow[:,line.line_id] .<= line.upstream_profile)
        @constraint(model, downstream_flow[:,line.line_id] .<= line.downstream_profile)
    end
    @expression(model, flow[i=1:NUM_TIMESTEPS,line=line_profiles.line_id], upstream_flow[i,line] .- downstream_flow[i,line])
    
    xprt = Dict{String, Vector{AffExpr}}()
    for zl in eachrow(zonal_limits)
        zone = zl.zone  
        xprt[zone] = Vector{AffExpr}(undef, NUM_TIMESTEPS)  
        for i in eachindex(xprt[zone])
            xprt[zone][i] = AffExpr(0.0)
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
        if !isempty(market[market.zone.==zl.zone,:])
            @constraint(model,
            sum(market[market.zone.==zl.zone,:].x .* market[market.zone.==zl.zone,:].profile, dims=1)[1] .>= xprt[zl.zone])
            @constraint(model,
            sum(market[market.zone.==zl.zone,:].x .* market[market.zone.==zl.zone,:].profile, dims=1)[1] .<= xprt[zl.zone])
        else
            @constraint(model,
            xprt[zl.zone] .== zeros(NUM_TIMESTEPS))
        end
    end
    @constraint(model,
        [b in eachrow(market)],
        b.s >= sum(b.profile .* (b.p .- b.price))
    )

    @constraint(model,
        [l in eachrow(line_profiles)],
        p[:,l.line_to] .- p[:,l.line_from] .<= upstream_profit[:,l.line_id]
    )

    @constraint(model,
        [l in eachrow(line_profiles)],
        p[:,l.line_from] .- p[:,l.line_to] .<= downstream_profit[:,l.line_id]
    )

    @constraint(
        model,
        duality_constraint,
        sum(market.s)
        + sum([sum(upstream_profit[:, line.line_id] .* line.upstream_profile) for line in eachrow(line_profiles)])
        + sum([sum(downstream_profit[:, line.line_id] .* line.downstream_profile) for line in eachrow(line_profiles)])
        == -sum(sum(market.x .* market.profile .* market.price))
    )

    @objective(model, Max, -sum(sum(market.x .* market.profile .* market.price))
        - 1e-7*(sum([sum(upstream_flow[:, line.line_id]) for line in eachrow(line_profiles)])
        + sum([sum(downstream_flow[:, line.line_id]) for line in eachrow(line_profiles)]))) # added penalty term to discourage abusing cycles TODO: add a post-solve model to simplify flows without this
    optimize!(model)
    @assert is_solved_and_feasible(model)
    market.accepted_quantity = [value(m.x).*m.profile for m in eachrow(market)]
    market.s = [value(m.s) for m in eachrow(market)]
    zonal_limits.zonal_price = [Array(value.(zl.p)) for zl in eachrow(zonal_limits)]
    zonal_limits.generation = [sum(
        reduce(
            hcat,
            [max.(aq, 0.0) for aq in market[market.zone .== zl.zone, :].accepted_quantity];
            init = zeros(24, 1)
        )',
        dims=1
    ) for zl in eachrow(zonal_limits)]
    zonal_limits.consumption = [sum(
        reduce(
            hcat, 
            [min.(aq, 0.0) for aq in market[market.zone .== zl.zone, :].accepted_quantity];
            init = zeros(24, 1)
        )',
        dims=1
    ) for zl in eachrow(zonal_limits)]

    line_profiles.flow = [Array(value.(flow[:,l.line_id])) for l in eachrow(line_profiles)]


    return (market[:, [:id, :zone, :accepted_quantity]],
            zonal_limits[:, [:zone, :zonal_price, :consumption, :generation]],
            line_profiles[:,[:line_from, :line_to, :flow]])
    
end

function recover_transfers(market_results::DataFrame)::DataFrame
end

# File output
function publish_results(market_results::DataFrame, transfers::DataFrame, filepath::String)
end

export load_lines
export load_bids
export Line
export BadFileError
export DuplicateLineError

end # module EasyDAM
