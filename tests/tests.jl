using EasyDAM
using Test

@test_throws EasyDAM.BadFileError load_lines("aaa")
@test_throws EasyDAM.BadFileError load_lines("./tests/resources/market_topology_ill_formed.yaml")
@test_throws EasyDAM.DuplicateLineError load_lines("./tests/resources/market_topology_with_duplicates.yaml")
@test length(load_lines("./tests/resources/market_topology.yaml"))==8

@test_throws EasyDAM.BadFileError load_bids("aaa")
@test_throws EasyDAM.BadFileError load_bids("./tests/resources/bids_ill_formed.yaml")
@test_throws EasyDAM.DuplicateBidError load_bids("./tests/resources/bids_with_duplicates.yaml")
@test length(load_bids("./tests/resources/bids.yaml"))==3