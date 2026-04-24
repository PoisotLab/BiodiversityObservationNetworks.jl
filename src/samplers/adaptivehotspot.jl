"""
    AdaptiveHotspot <: BONSampler

Sampling for hotspots. Takes an auxiliary variable (e.g. uncertainty) to target 
samples toward. 

Starts at the global maximum of the target/uncertainty surface. Subsequent points 
are chosen to maximize a trade-off between the target value and spatial diversity 
(measured via the determinant of a kernel matrix).

Requires a `CandidatePool` with features. The first feature row is used as the
target/uncertainty surface; subsequent rows are ignored.


# Fields
- `n::Int`: number of sites to select (default 50)
- `scale`: Matérn kernel range parameter ρ (default 1.0)
- `smoothness`: Matérn kernel smoothness ν (default 0.5, equivalent
  to an exponential kernel)

# References
- Andrade-Pacheco, R., et al. (2020) TODO
"""
@kwdef struct AdaptiveHotspot <: BONSampler
    n::Int = 50
    scale = 1.0
    smoothness = 0.5
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

function _sample(
    rng::AbstractRNG,
    sampler::AdaptiveHotspot, 
    cpool::CandidatePool
)
    N = cpool.n

    ρ, ν = sampler.scale, sampler.smoothness

    # Use first feature row as target/uncertainty surface
    uncertainty = cpool.features[1, :]

    # Precompute full Matérn kernel matrix
    coords = cpool.coordinates
    K_full = [_matérn(norm(coords[:, i] .- coords[:, j]), ρ, ν) for i in 1:N, j in 1:N]


    completed = zeros(Bool, N)

    # Start from the maximum of the uncertainty surface
    _, imax = findmax(uncertainty)
 
    selected = [imax]

    # Greedily add points maximizing sum of value and diversity score
    for i in 2:sampler.n
        best_score = -Inf
        best_s = NaN
        for j in eachindex(cpool.keys)
            if !completed[j]
                proposal_idxs = vcat(selected, j)
                K = K_full[proposal_idxs, proposal_idxs]
                score = uncertainty[j] + sqrt(log(i)) * _h(K)

                if score > best_score
                    best_score = score
                    best_s = j
                end
            end 
        end
        completed[best_s] = true
        push!(selected, best_s)
    end
    return selected
end


@testitem "We can use AdaptiveHotspot" begin
    bon = sample(AdaptiveHotspot(), rand(30,20))

    @test bon isa BiodiversityObservationNetwork
end
