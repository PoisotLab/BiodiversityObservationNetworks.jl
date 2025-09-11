"""
    AdaptiveHotspot

Adaptive hotspot sampling prioritizes high-uncertainty
regions while encouraging spatial diversity via a kernel-based criterion. This
implementation follows [Andrade-Pacheco2020FinHot](@cite): start
at the maximum of the uncertainty surface and iteratively add locations that
optimize a trade-off between local value and diversity with respect to already
chosen sites using a Matérn covariance kernel.

Arguments:
- `num_nodes`: number of sites to select
- `scale` (`ρ`): range parameter of the Matérn kernel
- `smoothness` (`ν`): smoothness parameter of the Matérn kernel
"""
@kwdef mutable struct AdaptiveHotspot{I <: Integer, F <: Real} <: BONSampler
    num_nodes::I = _DEFAULT_NUM_NODES
    scale::F = 1.
    smoothness::F = 0.5 # equivalent to exponential covariance  
end


# TODO:
# biased toward exterior points.
# for rasters, choose matern values from the (n,m) raster sized window around a particular point by reflecting the raster on the border, and compute distance matrix. i think average distance to all points should be equal then?


function _sample(
    sampler::AdaptiveHotspot, 
    uncertainty;
    inclusion = nothing,
)
    # Matérn kernel hyperparameters
    ρ, ν = sampler.scale, sampler.smoothness
    pool = getpool(uncertainty)
    features = getfeatures(uncertainty)

    done_idx = zeros(Bool, length(pool))

    # Pairwise distance matrix over candidate locations
    dist_mat = getdistancematrix(uncertainty)

    _M = _matérn.(dist_mat, ρ, ν)

    # Start from the maximum of the uncertainty surface
    imax = last(findmax(vec(features)))
 
    coords = [pool[imax]]
    selected_idxs = [imax]

    # Greedily add points maximizing sum of value and diversity score
    for i in 2:(sampler.num_nodes)
        best_score = -Inf
        best_s = NaN
        for (ci, cs) in enumerate(pool)
            if !done_idx[ci]
                proposal_idxs = vcat(selected_idxs, ci)
                K = _M[proposal_idxs, proposal_idxs]
                score = features[ci] + sqrt(log(i)) * _h(K)

                if score > best_score
                    best_score = score
                    best_s = ci
                end
            end 
        end
        done_idx[best_s] = true
        push!(selected_idxs, best_s)
        push!(coords, pool[best_s])
    end

    features = features[:, selected_idxs]
    return BiodiversityObservationNetwork(coords, features)
end

"""
    _matérn(d, ρ, ν)

Matérn covariance kernel evaluated at distance `d`, range `ρ`, smoothness `ν`.
Normalized so that `_matérn(0, ρ, ν) == 1`.
"""
function _matérn(d, ρ, ν)
    d == 0 && return 1
    return 1.0 * (2.0^(1.0 - ν)) / SpecialFunctions.gamma(ν) *
           (sqrt(2ν) * d / ρ)^ν *
           besselk(ν, sqrt(2ν) * d / ρ)
end

"""
    _h(K)

Entropy-like diversity score based on the log-determinant of the kernel matrix
`K`. Larger values encourage selecting points that are diverse with respect to
previously chosen sites.
"""
function _h(K)
    return (0.5 * log(2 * π * ℯ)^size(K,1)-1) * logabsdet(K)[1]
end
