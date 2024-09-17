"""
    SimpleRandom

Implements Simple Random spatial sampling (each location has equal probability of selection)
"""
Base.@kwdef struct SimpleRandom{I<:Integer} <: BONSampler
    numsites::I = 30
    function SimpleRandom(numsites::I) where{I<:Integer}
        srs = new{typeof(numsites)}(numsites)
        check_arguments(srs)
        return srs
    end
end

function check_arguments(srs::SRS) where {SRS <: SimpleRandom}
    check(TooFewSites, srs)
    return nothing
end

_default_pool(::SimpleRandom) = pool((50,50))

function _sample!(
    selections::S,
    candidates::C,
    sampler::SimpleRandom{I},
) where {S<:Sites,C<:Sites,I}
    _coords = Distributions.sample(candidates.coordinates, sampler.numsites; replace = false)
    for (i,c) in enumerate(_coords)
        selections[i] = c
    end
    return selections
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "SimpleRandom constructor with no arguments works" begin
    sr = SimpleRandom()
    @test typeof(sr) <: SimpleRandom
end

@testitem "SimpleRandom must have more than one point" begin
    @test_throws TooFewSites SimpleRandom(numsites = -1)
    @test_throws TooFewSites SimpleRandom(numsites = 0)
    @test_throws TooFewSites SimpleRandom(numsites = 1)
end

@testitem "SimpleRandom allows keyword arguments for number of points" begin
    N = 314
    srs = SimpleRandom(; numsites = N)
    @test srs.numsites == N
end
