"""
    BalancedAcceptance

A `BONSeeder` that uses Balanced-Acceptance Sampling (Van-dem-Bates et al. 2017
https://doi.org/10.1111/2041-210X.13003)
"""
Base.@kwdef struct BalancedAcceptance{I <: Integer} <: BONSeeder
    numpoints::I = 30
    dims::Tuple{I, I} = (50, 50)
    function BalancedAcceptance(numpoints, dims)
        bas = new{typeof(numpoints)}(numpoints, dims)
        check_arguments(bas)
        return bas
    end
end

function check_arguments(bas::BalancedAcceptance)
    check(TooFewSites, bas)
    max_num_sites = prod(bas.dims)
    return max_num_sites >= bas.numpoints || throw(
        TooManySites(
            "Number of sites to select $(bas.numpoints) is greater than number of possible sites $(max_num_sites)",
        ),
    )
end

function _generate!(
    coords::Vector{CartesianIndex},
    ba::BalancedAcceptance,
)
    seed = rand(Int32.(1e0:1e7), 2)
    x, y = ba.dims
    for (idx, ptct) in enumerate(eachindex(coords))
        i, j = haltonvalue(seed[1] + ptct, 2), haltonvalue(seed[2] + ptct, 3)
        coords[idx] = CartesianIndex(convert.(Int32, [ceil(x * i), ceil(y * j)])...)
    end
    return coords
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "BalancedAcceptance default constructor works" begin
    @test typeof(BalancedAcceptance()) <: BalancedAcceptance
end

@testitem "BalancedAcceptance requires positive number of sites" begin
    @test_throws TooFewSites BalancedAcceptance(numpoints = 1)
    @test_throws TooFewSites BalancedAcceptance(numpoints = 0)
    @test_throws TooFewSites BalancedAcceptance(numpoints = -1)
end

@testitem "BalancedAcceptance can't be run with too many sites" begin
    numpts, numcandidates = 26, 25
    @test numpts > numcandidates   # who watches the watchmen?
    dims = Int32.(floor.((sqrt(numcandidates), sqrt(numcandidates))))
    @test_throws TooManySites BalancedAcceptance(numpts, dims)
end

@testitem "BalancedAcceptance can generate points" begin
    bas = BalancedAcceptance()
    coords = seed(bas)

    @test typeof(coords) <: Vector{CartesianIndex}
    @test length(coords) == bas.numpoints
end

@testitem "BalancedAcceptance can generate a custom number of points as positional argument" begin
    numpts = 77
    sz = (50, 50)
    bas = BalancedAcceptance(numpts, sz)
    coords = seed(bas)
    @test numpts == length(coords)
end

@testitem "BalancedAcceptance can take number of points as keyword argument" begin
    N = 40
    bas = BalancedAcceptance(; numpoints = N)
    @test bas.numpoints == N
end
