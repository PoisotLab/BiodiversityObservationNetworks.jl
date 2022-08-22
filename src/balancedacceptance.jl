@kwdef mutable struct BalancedAcceptance{I <: Integer} <: BONSeeder
    numpoints::I = 50
    α = 1.0
end

function _validate(sampler::BalancedAcceptance)
    sampler.α < zero(sampler.α) || throw(ArgumentError("yada yada yada"))
    _checknumpoints(sampler)
    return nothing
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