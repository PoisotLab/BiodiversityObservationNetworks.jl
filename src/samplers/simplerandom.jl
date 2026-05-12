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
