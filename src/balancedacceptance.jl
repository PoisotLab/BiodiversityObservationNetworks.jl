"""
    BalancedAcceptance

A `BONSeeder` that uses Balanced-Acceptance Sampling (Van-dem-Bates et al. 2017
https://doi.org/10.1111/2041-210X.13003)
"""
Base.@kwdef mutable struct BalancedAcceptance{I <: Integer, F <: AbstractFloat} <: BONSeeder
    numpoints::I = 50
    α::F = 1.0
    function BalancedAcceptance(numpoints, α)
        if numpoints < one(numpoints)
            throw(
                ArgumentError(
                    "You cannot have a BalancedAcceptance with fewer than one point",
                ),
            )
        end
        if α <= zero(α)
            throw(
                ArgumentError(
                    "The value of α for BalancedAcceptance must be larger than 0",
                ),
            )
        end
        return new{typeof(numpoints), typeof(α)}(numpoints, α)
    end
end

function _generate!(
    coords::Vector{CartesianIndex},
    sampler::BalancedAcceptance,
    uncertainty::Matrix{T},
) where {T <: AbstractFloat}
    seed = rand(Int32.(1e0:1e7), 2)
    np, α = sampler.numpoints, sampler.α
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
            stduncert[i] = 1.
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
