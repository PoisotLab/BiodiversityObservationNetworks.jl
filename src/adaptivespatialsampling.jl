@kwdef struct AdaptiveSpatialSampling{I<:Integer} <: SpatialSampler
    numpoints::I = 50
end

function _generate!(ass::AdaptiveSpatialSampling, uncertainty::M) where {M<:AbstractMatrix}
    return coords
end

function H(threshold::T, uncertainty::Matrix{T}) where {T<:Number}
    p = mean(uncertainty .> threshold)
    q = 1.0 - p
    (isone(q) | iszero(q)) && return 0.0
    return -p * log2(p) - q * log2(q)
end

function matérn(d, ρ, ν, σ²)
    # ν = 0.5 to have the exponential version
    return σ² * (2.0^(1.0 - ν)) / gamma(ν) *
           (sqrt(2ν) * (d / ρ))^ν *
           bessely(ν, sqrt(2ν) * (d / ρ))
end

function h(d, ρ, ν, σ²)
    m₁ = length(d) # Number of cells in the current pool
    return (0.5 * log(2 * π * ℯ)^m₁) * sum(matérn(d, ρ, ν, σ²))
end

#=
using SpecialFunctions
using Statistics
using NeutralLandscapes

u = rand(DiamondSquare(), (20, 20))
pool = vcat(CartesianIndices(u)...)
s = eltype(pool)[]
imax = last(findmax([H(u[i], u) for i in pool]))
push!(s, popat!(pool, imax))
=#