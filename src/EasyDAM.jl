module EasyDAM
import YAML
import StructTypes
using DataFrames

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

function flatten_bids(bids::Vector{Bid})::DataFrame
    flattened_bids = DataFrame(id=UInt8[],price=Float64[],MAR=Float64[],zone=String[],profile=Vector{Float64}[])
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
function get_zonal_limits(lines::Vector{Line},zones::Vector{String})::DataFrame
    #zonal_limits = DataFrame(zone=String[],export_limit=Vector{Float64}[],import_limit=Vector{Float64}[])
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
    return DataFrame(zone=zones,export_limit=export_profiles,import_limit=import_profiles)
end

function build_market(bids::DataFrame,lines::DataFrame)::DataFrame
end

# Home of the magical gnomes
function solve_market(market::DataFrame)::DataFrame
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
