"""
    WeightedBalancedAcceptance

A `BONSeeder` that uses Balanced-Acceptance Sampling [@cite] combined with rejection sampling to create a set of sampling sites that is weighted toward values with higher uncertainty as a function of the parameter α.
"""
Base.@kwdef struct WeightedBalancedAcceptance{I <: Integer, F <: Real} <: BONSampler
    numsites::I = 30
    α::F = 1.0
    function WeightedBalancedAcceptance(numsites, α)
        wbas = new{typeof(numsites), typeof(α)}(numsites, α)
        check_arguments(wbas)
        return wbas
    end
end


function check_arguments(wbas::WeightedBalancedAcceptance)
    check(TooFewSites, wbas)
    wbas.α > 0 ||
        throw(
            ArgumentError("WeightedBalancedAcceptance requires α to be greater than 0 "),
        )
    return nothing
end

function _sample!(
    coords::S,
    candidates::C,
    sampler::WeightedBalancedAcceptance{I, T},
    weights::L
) where {S<:Sites,C<:Sites,I <: Integer, T <: AbstractFloat,L<:Layer}
    seed = rand(I.(1e0:1e7), 2)
    α = sampler.α
    x, y = size(weights)

    nonnan_indices = findall(!isnan, weights)
    stduncert = similar(weights)

    uncert_values = weights[nonnan_indices]
    stduncert_values = similar(uncert_values)
    zfit = nothing
    if var(uncert_values) > 0
        zfit = StatsBase.fit(ZScoreTransform, uncert_values)
        stduncert_values = StatsBase.transform(zfit, uncert_values)
    end

    nonnan_counter = 1
    for i in eachindex(weights)
        if isnan(weights[i])
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

    return coords
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
    @test_throws TooFewSites WeightedBalancedAcceptance(numsites = 0)
    @test_throws TooFewSites WeightedBalancedAcceptance(numsites = 1)
end

@testitem "WeightedBalancedAcceptance can take bias parameter α as keyword argument" begin
    α = 3.14159
    wbas = WeightedBalancedAcceptance(; α = α)
    @test wbas.α == α
end

@testitem "WeightedBalancedAcceptance can take number of points as keyword argument" begin
    N = 40
    wbas = WeightedBalancedAcceptance(; numsites = N)
    @test wbas.numsites == N
end
