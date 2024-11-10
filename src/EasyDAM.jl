module EasyDAM
import YAML
import StructTypes

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
    if !(:lines in key(lines_dict))
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
    if !(:bids in key(bids_dict))
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

export load_lines
export load_bids
export Line
export BadFileError
export DuplicateLineError

end # module EasyDAM
