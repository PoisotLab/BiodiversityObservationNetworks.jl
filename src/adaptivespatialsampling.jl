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
           besselk(ν, sqrt(2ν) * (d / ρ))
end

function h(d, ρ, ν, σ²)
    K = [matérn(i, ρ, ν, σ²) for i in d]
    return (0.5 * log(2 * π * ℯ)^length(d)) * sum(filter(!isnan, K))
end

#=
using SpecialFunctions
using Statistics
using Plots
using NeutralLandscapes
using Plots

function D(a1::T, a2::T) where {T <: CartesianIndex{2}}
    x1, y1 = first(a1.I), first(a2.I)
    x2, y2 = last(a1.I), last(a2.I)
    return sqrt((x1-x2)^2.0+(y1-y2)^2.0)
end

u = rand(PerlinNoise((4,4)), (60, 60))
heatmap(u, c=:viridis)

pool = vcat(CartesianIndices(u)...)
s = eltype(pool)[]

imax = last(findmax([u[i] for i in pool]))
push!(s, popat!(pool, imax))
scatter!([reverse(x.I) for x in s], lab="", c=:white)

for t in 1:10
    candidates_s = [push!(copy(s), p) for p in pool]
    st = zeros(Float32, length(candidates_s))
    for (ci, cs) in enumerate(candidates_s)
        d = reduce(vcat, [[D(cs[i], cs[j]) for j in (i+1):length(cs)] for i in 1:(length(cs)-1)])
        st[ci] = u[last(cs)] + sqrt(log(t)) * h(d, 1.0, 0.5, var(u[cs]))
    end
    push!(s, popat!(pool, last(findmax(st))))
end
scatter!([reverse(x.I) for x in s], lab="", c=:white)
=#