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

function D(a1::T, a2::T) where {T <: CartesianIndex{2}}
    x1, y1 = a1.I
    x2, y2 = a2.I
    return sqrt((x1-x2)^2.0+(y1-y2)^2.0)
end

#=
using SpecialFunctions
using Statistics
using Plots
using NeutralLandscapes

u = rand(PerlinNoise((4,4)), (60, 60))
heatmap(u, c=:viridis, cbar=false, frame=:none, aspectratio=1, dpi=500)

pool = vcat(CartesianIndices(u)...)
s = Vector{eltype(pool)}(undef, 250)

imax = last(findmax([u[i] for i in pool]))
s[1] = popat!(pool, imax)

@time for t in 2:length(s)
    best_score = 0.0
    best_s = 1
    for (ci, cs) in enumerate(pool)
        s[t] = cs
        d = reduce(vcat, [[D(s[i], s[j]) for j in (i+1):t] for i in 1:(t-1)])
        score = u[cs] + sqrt(log(t)) * h(d, 1.0, 0.5, var(u[s[1:t]]))
        if score > best_score
            best_score = score
            best_s = ci
        end
    end
    s[t] = popat!(pool, best_s)
end

scatter!([reverse(x.I) for x in s], lab="", c=:white, cbar=false)

savefig(joinpath(homedir(), "lol.png"))
=#