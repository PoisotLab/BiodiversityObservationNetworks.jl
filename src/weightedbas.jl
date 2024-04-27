"""
    WeightedBalancedAcceptance

A `BONSeeder` that uses Balanced-Acceptance Sampling [@cite] combined with rejection sampling to create a set of sampling sites that is weighted toward values with higher uncertainty as a function of the parameter α.
"""
Base.@kwdef struct WeightedBalancedAcceptance{I <: Integer, F <: Real} <: BONSeeder
    numpoints::I = 3
    α::F = 1.0
    function WeightedBalancedAcceptance(numpoints, α)
        wbas = new{typeof(numpoints), typeof(α)}(numpoints, α)
        _check_arguments(wbas)
        return wbas
    end
end

function _generate!(
    coords::Vector{CartesianIndex},
    sampler::WeightedBalancedAcceptance,
    uncertainty::Matrix{T},
) where {T <: AbstractFloat}
    seed = rand(Int32.(1e0:1e7), 2)
    α = sampler.α
    x, y = size(uncertainty)

    nonnan_indices = findall(!isnan, uncertainty)
    stduncert = similar(uncertainty)

    uncert_values = uncertainty[nonnan_indices]
    stduncert_values = similar(uncert_values)
    zfit = nothing
    if var(uncert_values) > 0
        zfit = StatsBase.fit(ZScoreTransform, uncert_values)
        stduncert_values = StatsBase.transform(zfit, uncert_values)
    end

    nonnan_counter = 1
    for i in eachindex(uncertainty)
        if isnan(uncertainty[i])
            stduncert[i] = NaN
        elseif !isnothing(zfit)
            stduncert[i] = stduncert_values[nonnan_counter]
            nonnan_counter += 1
        else
            stduncert[i] = 1.0
        end
    end

    reluncert = broadcast(x -> isnan(x) ? NaN : exp(α * x) / (1 + exp(α * x)), stduncert)
    ptct = 1
    addedpts = 1
    while addedpts <= length(coords)
        i, j = haltonvalue(seed[1] + ptct, 2), haltonvalue(seed[2] + ptct, 3)
        candcoord = CartesianIndex(convert.(Int32, [ceil(x * i), ceil(y * j)])...)
        prob = reluncert[candcoord]
        if !isnan(prob) && rand() < prob
            coords[addedpts] = candcoord
            addedpts += 1
        end
        ptct += 1
    end

    return (coords, uncertainty)
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "WeightedBalancedAcceptance default constructor works" begin
    @test typeof(WeightedBalancedAcceptance()) <: WeightedBalancedAcceptance
end

@testitem "WeightedBalancedAcceptance requires positive number of sites" begin
    α = 1.0
    @test_throws TooFewSites WeightedBalancedAcceptance(0, α)
    @test_throws TooFewSites WeightedBalancedAcceptance(1, α)
end

@testitem "WeightedBalancedAcceptance can't be run with too many sites" begin
    α = 1.0
    numpts, numcandidates = 26, 25
    @test numpts > numcandidates   # who watches the watchmen?
    wbas = WeightedBalancedAcceptance(numpts, α)
    dims = Int32.(floor.([sqrt(numcandidates), sqrt(numcandidates)]))
    uncert = rand(dims...)

    @test_throws TooManySites seed(wbas, uncert)
end

@testitem "WeightedBalancedAcceptance can generate points" begin
    wbas = WeightedBalancedAcceptance()
    sz = (50, 50)
    coords = seed(wbas, rand(sz...)) |> first

    @test typeof(coords) <: Vector{CartesianIndex}
    @test length(coords) == wbas.numpoints
end

@testitem "WeightedBalancedAcceptance can generate a custom number of points" begin
    numpts = 77
    α = 1.0
    wbas = WeightedBalancedAcceptance(numpts, α)
    sz = (50, 50)
    coords = seed(wbas, rand(sz...)) |> first
    @test numpts == length(coords)
end

@testitem "BalancedAcceptance can take bias parameter α as keyword argument" begin
    α = 3.14159
    wbas = WeightedBalancedAcceptance(; α = α)
    @test wbas.α == α
end

@testitem "BalancedAcceptance can take number of points as keyword argument" begin
    N = 40
    wbas = WeightedBalancedAcceptance(; numpoints = N)
    @test wbas.numpoints == N
end
