"""
    SimpleRandom

Implements Simple Random spatial sampling (each location has equal probability of selection)
"""
Base.@kwdef struct SimpleRandom{I <: Integer} <: BONSeeder
    numpoints::I = 30
    dims::Tuple{I, I} = (50, 50)
    function SimpleRandom(numpoints, dims)
        srs = new{typeof(numpoints)}(numpoints, dims)
        check_arguments(srs)
        return srs
    end
end

function check_arguments(srs::SimpleRandom)
    check(TooFewSites, srs)
    max_num_sites = prod(srs.dims)
    return max_num_sites >= srs.numpoints || throw(
        TooManySites(
            "Number of sites to select $(srs.numpoints) is greater than number of possible sites $(max_num_sites)",
        ),
    )
end

function _generate!(
    coords::Vector{CartesianIndex},
    sampler::SimpleRandom,
)
    pool = CartesianIndices(sampler.dims)
    coords .= sample(pool, sampler.numpoints; replace = false)
    return coords
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
    @test_throws TooFewSites SimpleRandom(numpoints = -1)
    @test_throws TooFewSites SimpleRandom(numpoints = 0)
    @test_throws TooFewSites SimpleRandom(numpoints = 1)
end

@testitem "SimpleRandom allows keyword arguments for number of points" begin
    N = 314
    srs = SimpleRandom(; numpoints = N)
    @test srs.numpoints == N
end

@testitem "SimpleRandom throws exception if there are more sites than candidates" begin
    numpts, numcandidates = 26, 25
    dims = Int.(floor.((sqrt(numcandidates), sqrt(numcandidates))))
    srs = @test prod(dims) < numpts
    @test_throws TooManySites SimpleRandom(; numpoints = numpts, dims = dims)
end
