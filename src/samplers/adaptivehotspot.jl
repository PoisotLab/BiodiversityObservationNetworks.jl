"""
    AdaptiveHotspot

@Andrade-Pacheco2020
"""
Base.@kwdef mutable struct AdaptiveHotspot{I <: Integer} <: BONSampler
    number_of_nodes::I = 50
end

_valid_geometries(::AdaptiveHotspot) = (Raster)


function _sample(sampler::AdaptiveHotspot, uncertainty::Raster, bon::BiodiversityObservationNetwork)
    
end


# TODO: turn this pool call into separate dispatch for Raster's, and then have
# the internal _adaptive_hotspot work on BONs, much faster because
# pairwise distance matrix takes forever for full rasters.
function _sample(sampler::AdaptiveHotspot, uncertainty::Raster)
    pool = nonempty(uncertainty)
    d = zeros(Float64, Int((sampler.number_of_nodes * (sampler.number_of_nodes - 1)) / 2))
    imax = last(findmax([uncertainty[i] for i in pool]))
 
    coords = Array{CartesianIndex,1}(undef, sampler.number_of_nodes)
    coords[begin] = popat!(pool, imax)

    best_score = 0.0
    best_s = 1

    for i in 2:(sampler.number_of_nodes)
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

    Es, Ns = SDT.eastings(uncertainty), SDT.northings(uncertainty)

    return BiodiversityObservationNetwork([Node(Es[c[2]], Ns[c[1]]) for c in coords])
end


function _matérn(d, ρ, ν)
    return 1.0 * (2.0^(1.0 - ν)) / SpecialFunctions.gamma(ν) *
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

