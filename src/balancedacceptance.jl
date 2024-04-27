"""
    BalancedAcceptance

A `BONSeeder` that uses Balanced-Acceptance Sampling (Van-dem-Bates et al. 2017
https://doi.org/10.1111/2041-210X.13003)
"""
Base.@kwdef struct BalancedAcceptance{I <: Integer} <: BONSeeder
    numpoints::I = 30
    function BalancedAcceptance(numpoints)
        bas = new{typeof(numpoints)}(numpoints)
        _check_arguments(bas)
        return bas
    end
end

function _generate!(
    coords::Vector{CartesianIndex},
    ::BalancedAcceptance,
    uncertainty::Matrix{T},
) where {T <: Real}
    seed = rand(Int32.(1e0:1e7), 2)
    x, y = size(uncertainty)
    for (idx, ptct) in enumerate(eachindex(coords))
        i, j = haltonvalue(seed[1] + ptct, 2), haltonvalue(seed[2] + ptct, 3)
        coords[idx] = CartesianIndex(convert.(Int32, [ceil(x * i), ceil(y * j)])...)
    end
    return (coords, uncertainty)
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
    @test_throws TooFewSites BalancedAcceptance(0)
    @test_throws TooFewSites BalancedAcceptance(1)
end

@testitem "BalancedAcceptance can't be run with too many sites" begin
    numpts, numcandidates = 26, 25
    @test numpts > numcandidates   # who watches the watchmen?
    bas = BalancedAcceptance(numpts)
    dims = Int32.(floor.([sqrt(numcandidates), sqrt(numcandidates)]))
    uncert = rand(dims...)

    @test_throws TooManySites seed(bas, uncert)
end

@testitem "BalancedAcceptance can generate points" begin
    bas = BalancedAcceptance()
    sz = (50, 50)
    coords = seed(bas, rand(sz...)) |> first

    @test typeof(coords) <: Vector{CartesianIndex}
    @test length(coords) == bas.numpoints
end

@testitem "BalancedAcceptance can generate a custom number of points" begin
    numpts = 77
    bas = BalancedAcceptance(numpts)
    sz = (50, 50)
    coords = seed(bas, rand(sz...)) |> first
    @test numpts == length(coords)
end

@testitem "BalancedAcceptance can take number of points as keyword argument" begin
    N = 40
    bas = BalancedAcceptance(; numpoints = N)
    @test bas.numpoints == N
end
