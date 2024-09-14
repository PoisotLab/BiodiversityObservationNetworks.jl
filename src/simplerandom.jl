"""
    SimpleRandom

Implements Simple Random spatial sampling (each location has equal probability of selection)
"""
Base.@kwdef struct SimpleRandom{I<:Integer,S<:Sites} <: BONSampler
    numsites::I = 30
    pool::S = Sites(vec(collect(CartesianIndices((1:50, 1:50)))))
    function SimpleRandom(numsites::I, pool::J) where{I<:Integer,J<:Sites}
        srs = new{typeof(numsites), typeof(pool)}(numsites, pool)
        check_arguments(srs)
        return srs
    end
end
maxsites(srs::SimpleRandom) = length(pool(srs))

function SimpleRandom(layer::Layer, numsites = 50)
    candidates = pool(layer)
    srs = SimpleRandom(numsites, candidates)
    check_arguments(srs)
    return srs
end


function check_arguments(srs::SRS) where {SRS <: SimpleRandom}
    check(TooFewSites, srs)
    check(TooManySites, srs)
    return nothing
end

function _sample!(
    sites::Sites{T},
    sampler::SimpleRandom{I},
) where {T,I}
    candidates = coordinates(pool(sampler))
    _coords = Distributions.sample(candidates, sampler.numsites; replace = false)
    for (i,c) in enumerate(_coords)
        sites[i] = c
    end
    return sites
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

@testitem "SimpleRandom throws exception if there are more sites than candidates" begin
    numpts, numcandidates = 26, 25
    dims = Int.(floor.((sqrt(numcandidates), sqrt(numcandidates))))
    srs = @test prod(dims) < numpts
    @test_throws TooManySites SimpleRandom(; numsites = numpts, dims = dims)
end
