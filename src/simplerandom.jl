"""
    SimpleRandom

Implements Simple Random spatial sampling (each location has equal probability of selection)
"""
Base.@kwdef struct SimpleRandom{I <: Integer} <: BONSeeder
    numpoints::I = 50
    function SimpleRandom(numpoints)
        srs = new{typeof(numpoints)}(numpoints)
        _check_arguments(srs)
        return srs
    end
end

function _generate!(
    coords::Vector{CartesianIndex},
    sampler::SimpleRandom,
    uncertainty::Matrix{T},
) where {T <: AbstractFloat}
    pool = CartesianIndices(uncertainty)

    coords .= sample(pool, sampler.numpoints; replace = false)
    return (coords, uncertainty)
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
    @test_throws TooFewSites SimpleRandom(-1)
    @test_throws TooFewSites SimpleRandom(0)
    @test_throws TooFewSites SimpleRandom(1)
end

@testitem "SimpleRandom allows keyword arguments for number of points" begin
    N = 314
    srs = SimpleRandom(; numpoints = N)
    @test srs.numpoints == N
end

@testitem "SimpleRandom throws exception if there are more sites than candidates" begin
    numpts, numcandidates = 26, 25
    srs = SimpleRandom(; numpoints = numpts)
    dims = Int32.(floor.([sqrt(numcandidates), sqrt(numcandidates)]))
    uncert = rand(dims...)
    @test prod(dims) < numpts
    @test_throws TooManySites seed(srs, uncert)
end
