@kwdef struct AdaptiveSpatialSampling{I<:Integer} <: SpatialSampler
    numpoints::I = 50
end

function _generate!(ass::AdaptiveSpatialSampling, uncertainty::M) where {M<:AbstractMatrix}
    return coords
end

function H(threshold::T, uncertainty::Matrix{T}) where {T<:Number}
    p = mean(uncertainty .> threshold)
    q = 1.0 - p
    return -p * log2(p) - q * log2(q)
end

function matérn(d, ρ, ν, σ²)
    return σ² * (2.0^(1.0 - ν)) / gamma(ν) *
           (sqrt(2ν) * (d / ρ))^ν *
           bessely(ν, sqrt(2ν) * (d / ρ))
end
