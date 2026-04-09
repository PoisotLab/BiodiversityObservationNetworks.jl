"""
    SimpleRandom <: BONSampler

Simple random sampling (or weighted random sampling when inclusion weights
are non-uniform). Selects `n` candidates without replacement.
"""
@kwdef struct SimpleRandom <: BONSampler
    n::Int = 50
end

supports_inclusion(::SimpleRandom) = true
guarantees_exact_n(::SimpleRandom) = true

function _sample(rng::AbstractRNG, sampler::SimpleRandom, cpool::CandidatePool)
    weights = StatsBase.Weights(cpool.inclusion)
    selected = StatsBase.sample(rng, 1:cpool.n, weights, sampler.n; replace=false)
    return selected
end

@testitem "SimpleRandom works" begin
    result = sample(SimpleRandom(10), rand(20, 20))
    @test length(result) == 10
    @test length(unique(result.sites)) == 10
end

@testitem "SimpleRandom respects a mask" begin
    mask = falses(20, 20)
    mask[10:20, 10:20] .= true
    result = sample(SimpleRandom(10), rand(20, 20); mask)
    @test all(s -> s[1] >= 10 && s[2] >= 10, result.sites)
end 

#=

"""
    SimpleRandom <: BONSampler

Implements simple random sampling (SRS) and weighted random sampling.

# Fields
- `num_nodes::Int`: The number of sites to select.

# Description
Selects `num_nodes` locations from the domain uniformly at random without replacement. 
If an `inclusion` probability surface is provided, it performs weighted random sampling 
without replacement, where the weight of each candidate cell is proportional to its 
inclusion probability.

# Notes
While computationally efficient, Simple Random Sampling does not guarantee spatial 
balance and may result in clustering of sampling points.
"""
@kwdef struct SimpleRandom <: BONSampler
    num_nodes = _DEFAULT_NUM_NODES
end

function _sample(
    sampler::SimpleRandom,
    domain;
    inclusion = nothing
)
    pool = getpool(domain)

    nodes = isnothing(inclusion) ? SB.sample(pool, sampler.num_nodes, replace=false) : SB.wsample(pool, inclusion[pool], sampler.num_nodes, replace=false)
    return nodes, nothing
end 

# ========================================================================
# Tests
# ========================================================================

@testitem "We can use SimpleRandom with a RasterDomain" begin
    bon = sample(SimpleRandom(), rand(30,20))

    @test bon isa BiodiversityObservationNetwork
    @test first(bon) isa CartesianIndex
end


@testitem "We can use SimpleRandom with a RasterDomain and custom inclusion probabilities" begin
    inclusion = rand(30, 20)
    bon = sample(SimpleRandom(100), rand(30,20), inclusion=inclusion)

    @test bon isa BiodiversityObservationNetwork
    @test first(bon) isa CartesianIndex
end


@testitem "We can use SimpleRandom with a BON" begin
    candidate_bon = sample(SimpleRandom(100), rand(30,20))

    bon = sample(SimpleRandom(50), candidate_bon)

    @test bon isa BiodiversityObservationNetwork
    @test first(bon) isa CartesianIndex
end


@testitem "We can use SimpleRandom with a BON and custom inclusion probabilities" begin
    candidate_bon = sample(SimpleRandom(100), rand(30,20))
    inclusion = rand(100)
    bon = sample(SimpleRandom(50), candidate_bon, inclusion=inclusion)

    @test bon isa BiodiversityObservationNetwork
    @test first(bon) isa CartesianIndex
end

=#
