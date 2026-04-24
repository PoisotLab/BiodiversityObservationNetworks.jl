"""
    SpatiallyCorrelatedPoisson <: BONSampler

Implements Spatially Correlated Poisson Sampling (SCPS).

# Fields
- `n::Int`: expected number of sites to select (default 50)

# Description
Iterates through units, selecting them based on inclusion probabilities, and 
dynamically adjusting the probabilities of neighboring units to maintain spatial balance.

# References
- TODO set up refs
"""
@kwdef struct SpatiallyCorrelatedPoisson <: BONSampler
    n::Int = 50
end

supports_inclusion(::SpatiallyCorrelatedPoisson) = true

# Eq. (4) in Grafström (2012) bounds the maximum weight like this
_scps_max_weight(πᵢ, πⱼ) = min(πᵢ / (1.0 - πⱼ), (1.0 - πᵢ) / πⱼ)

function _sample(rng::AbstractRNG, sampler::SpatiallyCorrelatedPoisson, cpool::CandidatePool)
    N = cpool.n
    tree = _build_kdtree(cpool)

    # Scale inclusion to sum to n
    inclusion_probs = cpool.inclusion .* sampler.n

    inclusion_flags = zeros(Bool, N)
    prev_inclusion_probs = deepcopy(inclusion_probs)
    neighbor_weights = zeros(N)

    for j in 1:N
        # Step 1: Include unit j with probability inclusion_probs[j]
        inclusion_flags[j] = rand(rng) < inclusion_probs[j]
        flag = inclusion_flags[j]

        # Step 2: Reset weights and store current inclusion_probs
        prev_inclusion_probs .= inclusion_probs
        neighbor_weights .= 0.0
        πⱼ_prev = prev_inclusion_probs[j]

        # Step 3: Compute max allowable weights for remaining indices
        neighbors = _neighbor_order(tree, cpool.coordinates, j)
        for i in neighbors
            if i > j
                πᵢ = prev_inclusion_probs[i]
                neighbor_weights[i] = _scps_max_weight(πᵢ, πⱼ_prev)
            end 
        end

        # Step 4: Normalize weights so they sum to one
        w_sum = sum(neighbor_weights)
        if w_sum > 0
            neighbor_weights ./= w_sum
        end

        # Step 5: Update inclusion probabilities for unprocessed units (i.e. with index > j)
        for i in neighbors
            if i > j 
                inclusion_probs[i] = prev_inclusion_probs[i] - (flag - πⱼ_prev) * neighbor_weights[i]
            end
        end
    end

    return findall(inclusion_flags)
end

@testitem "SCPS runs" begin
    result = sample(SpatiallyCorrelatedPoisson(10), rand(20, 20))
    @test length(result) > 0
    @test all(s -> s isa CartesianIndex{2}, result.sites)
end


@testitem "SCPS works with a mask" begin
    mask = falses(20, 20)
    mask[1:10, 1:10] .= true
    result = sample(SpatiallyCorrelatedPoisson(10), rand(20, 20); mask)
    @test all(s -> s[1] <= 10 && s[2] <= 10, result.sites)
end

