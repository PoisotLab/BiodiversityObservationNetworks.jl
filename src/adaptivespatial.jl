"""
    AdaptiveSpatial

...

**numpoints**, an Integer (def. 50), specifying the number of points to use.

**α**, an AbstractFloat (def. 1.0), specifying ...
"""
Base.@kwdef mutable struct AdaptiveSpatial{T <: Integer} <: BONRefiner
    numpoints::T = 50
    function AdaptiveSpatial(numpoints)
        if numpoints < one(numpoints)
            throw(
                ArgumentError(
                    "You cannot have an AdaptiveSpatial with fewer than one point",
                ),
            )
        end
        return new{typeof(numpoints)}(numpoints)
    end
end

function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::AdaptiveSpatial,
    uncertainty::Matrix{T},
) where {T <: AbstractFloat}

    # Distance matrix (inlined)
    d = zeros(Float64, Int((sampler.numpoints * (sampler.numpoints - 1)) / 2))

    # Start with the point with maximum entropy
    imax = last(findmax([uncertainty[i] for i in pool]))
    # Add it to the stack
    coords[1] = popat!(pool, imax)

    best_score = 0.0
    best_s = 1

    for i in 2:(sampler.numpoints)
        for (ci, cs) in enumerate(pool)
            coords[i] = cs
            # Distance update
            start_from = Int((i - 1) * (i - 2) / 2) + 1
            end_at = start_from + Int(i - 2)
            d_positions = start_from:end_at
            for ti in 1:(i - 1)
                d[d_positions[ti]] = _D(cs, coords[ti])
            end
            # Get the score
            score = uncertainty[cs] + sqrt(log(i)) * _h(d[1:end_at], 1.0, 0.5)
            if score > best_score
                best_score = score
                best_s = ci
            end
        end
        coords[i] = popat!(pool, best_s)
    end
    return (coords, uncertainty)
end

function _H(threshold::T, uncertainty::Matrix{T}) where {T <: AbstractFloat}
    p = mean(uncertainty .> threshold)
    q = 1.0 - p
    (isone(q) | iszero(q)) && return 0.0
    return -p * log2(p) - q * log2(q)
end

function _matérn(d, ρ, ν)
    # This is the version from the supp mat
    # ν = 0.5 to have the exponential version
    return 1.0 * (2.0^(1.0 - ν)) / gamma(ν) *
           (sqrt(2ν) * d / ρ)^ν *
           besselk(ν, sqrt(2ν) * d / ρ)
end

function _h(d, ρ, ν)
    K = [_matérn(i, ρ, ν) for i in d]
    return (0.5 * log(2 * π * ℯ)^length(d)) * sum(K)
end

function _D(a1::T, a2::T) where {T <: CartesianIndex{2}}
    x1, y1 = a1.I
    x2, y2 = a2.I
    return sqrt((x1 - x2)^2.0 + (y1 - y2)^2.0)
end
