"""
    BalancedAcceptance

...

**numpoints**, an Integer (def. 50), specifying the number of points to use.

**α**, an AbstractFloat (def. 1.0), specifying ...
"""
Base.@kwdef mutable struct BalancedAcceptance{I <: Integer, F <: AbstractFloat} <: BONSeeder
    numpoints::I = 50
    α::F = 1.0
    function BalancedAcceptance(numpoints, α)
        if numpoints < one(numpoints)
            throw(ArgumentError("You cannot have a BalancedAcceptance with fewer than one point"))
        end
        if α <= zero(α)
            throw(ArgumentError("The value of α for BalancedAcceptance must be larger than 0"))
        end
        return new{typeof(numpoints), typeof(α)}(numpoints, α)
    end
end

function _generate!(coords::Vector{CartesianIndex}, sampler::BalancedAcceptance, uncertainty::Matrix{T}) where {T<:AbstractFloat}
    seed = Int32.([floor(10^7 * rand()), floor(10^7 * rand())])
    np, α = sampler.numpoints, sampler.α
    x, y = size(uncertainty)

    stduncert = StatsBase.transform(StatsBase.fit(ZScoreTransform, uncertainty, dims=2), uncertainty)
    reluncert = broadcast(x -> exp(α * x) / (1 + exp(α * x)), stduncert)

    ptct = 0
    while length(coords) < np
        i, j = haltonvalue(seed[1] + ptct, 2), haltonvalue(seed[2] + ptct, 3)
        candcoord = CartesianIndex(convert.(Int32, [ceil(x * i), ceil(y * j)])...)
        prob = reluncert[candcoord]
        if  rand() < prob 
            coords[ptct] = candcoord
            ptct += 1
        end 
    end

    return coords
end