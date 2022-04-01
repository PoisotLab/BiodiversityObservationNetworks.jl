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

function matérn(d, ρ, ν)
    # This is the version from the supp mat
    # ν = 0.5 to have the exponential version
    return 1.0 * (2.0^(1.0 - ν)) / gamma(ν) *
           (sqrt(2ν) * d / ρ)^ν *
           besselk(ν, sqrt(2ν) * d / ρ)
end

function h(d, ρ, ν)
    K = [matérn(i, ρ, ν) for i in d]
    return (0.5 * log(2 * π * ℯ)^length(d)) * sum(K)
end

function D(a1::T, a2::T) where {T <: CartesianIndex{2}}
    x1, y1 = a1.I
    x2, y2 = a2.I
    return sqrt((x1-x2)^2.0+(y1-y2)^2.0)
end

using SpecialFunctions
using Statistics
using Plots
using NeutralLandscapes
using StatsBase

u = rand(DiamondSquare(), (80, 80))
heatmap(u, c=:viridis, cbar=false, frame=:none, aspectratio=1, dpi=500)

pool = vcat(CartesianIndices(u)...)
pool = sample(pool, 500, replace=false)

scatter!([reverse(x.I) for x in pool], lab="", c=:black, ms=1)

steps = 75
s = Vector{eltype(pool)}(undef, steps)
d = zeros(Float64, Int((steps*(steps-1))/2))

imax = last(findmax([u[i] for i in pool]))
s[1] = popat!(pool, imax)

@time for t in 2:length(s)
    best_score = 0.0
    best_s = 1
    for (ci, cs) in enumerate(pool)
        s[t] = cs
        # Distance update
        start_from = Int((t-1)*(t-2)/2)+1
        end_at = start_from+Int(t-2)
        d_positions = start_from:end_at
        for ti in 1:(t-1)
            d[d_positions[ti]] = D(cs, s[ti])
        end
        # Get the score
        score = u[cs] + sqrt(log(t)) * h(d[1:end_at], 1.0, 0.5)
        if score > best_score
            best_score = score
            best_s = ci
        end
    end
    s[t] = popat!(pool, best_s)
end

scatter!([reverse(x.I) for x in s], lab="", c=:white, cbar=false)
