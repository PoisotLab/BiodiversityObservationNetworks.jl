"""
    SpatiallyCorrelatedPoisson

Implements the Spatially Correlated Poisson Sampling algorithm from Grafström (2012), Journal of Statistical Planning and Inference. doi:10.1016/j.jspi.2011.07.003
    
Iterate through all potential nodes i. Select it w.p. πᵢ, and adjust all πⱼ for all subsequent nodes j to maintain expected sample size.

I think the issue w/ samples occasionally being off-by-one is the last sample is rejected and the remaining inclusion probability has nowhere left to go. Same for Pivotal.
"""
@kwdef mutable struct SpatiallyCorrelatedPoisson <: BONSampler
    num_nodes = _DEFAULT_NUM_NODES
end

# Eq. (4) in Grafström (2012) bounds w_ij as the following
_max_weight(πᵢ, πⱼ) = min(πᵢ / (1 - πⱼ), (1 - πᵢ) / (πⱼ))

function _store!(rd1::RasterDomain, rd2::RasterDomain)
    rd1.data .= rd2.data
end 
_store!(v1, v2) = v1 .= v2

"""
    _sample(sampler::SpatiallyCorrelatedPoisson, domain, inclusion_probs)

Perform Spatially Correlated Poisson sampling following Grafström (2012).

# Arguments
- `sampler::SpatiallyCorrelatedPoisson`: The sampling strategy object.
- `domain`: The spatial domain (supports `size(domain)` and distance calculation).
- `inclusion_probs`: Initial inclusion probabilities for each unit.
        Defaults to uniform: `(number_of_nodes / N)` for all units.

# Returns
- A [`BiodiversityObservationNetwork.jl`](@ref) containing units included in the sample.

# Example Usage
`sample(SpatiallyCorrelatedPoisson(), zeros(30,20))`

"""
function _sample(
    sampler::SpatiallyCorrelatedPoisson,
    domain;
    inclusion=nothing,
)
    pool = getpool(domain)


    # If no custom inclusion probabilities are provided, use uniform inclusion
    inclusion_probs = isnothing(inclusion) ? get_uniform_inclusion(sampler, domain) : inclusion

    # Boolean inclusion indicator: true if unit is selected
    inclusion_flags = zeros(Bool, length(pool))

    # Copy of inclusion probabilities from previous iteration (j-1)
    prev_inclusion_probs = deepcopy(inclusion_probs)

    # Weight vector used for redistributing inclusion probabilities
    neighbor_weights = zeros(length(pool))

    # TODO: have option to truncate weights at k-th nearest neighbor
    # Neighbor order for each unit, excluding itself
    # neighbor_order[i] => [index of unit closest to i, index of unit 2nd closest to i, ...]
    _, neighbor_order = getnearestneighbors(domain)


    # Main loop: iterate over each unit j
    for j in 1:length(pool)
        # Step 1: Include unit j with probability inclusion_probs[j]

        πⱼ = inclusion_probs[pool[j]]
        inclusion_flags[j] = rand() < πⱼ
        inclusion_flag = inclusion_flags[j]

        # Step 2: Reset weights and store current inclusion_probs
        neighbor_weights .= 0
        _store!(prev_inclusion_probs, inclusion_probs)

        πⱼ_prev = prev_inclusion_probs[pool[j]]

        # Step 3: Compute max allowable weights (Eq. 4 constraint)
        for i in neighbor_order[j]
            πᵢ_prev = prev_inclusion_probs[pool[i]]

            # Ignore neighbors that have already has inclusion chosen, i.e. only consider units i > j 
            if i > j
                neighbor_weights[i] = _max_weight(πᵢ_prev, πⱼ_prev)
            end
        end

        # Step 4: Normalize weights so they sum to one
        neighbor_weights ./= sum(neighbor_weights)

        # Step 5: Update inclusion_probs for units i > j
        for i in neighbor_order[j]
            πᵢ_prev = prev_inclusion_probs[pool[i]]

            if i > j
                # Apply update rule
                inclusion_probs[i] = πᵢ_prev - (inclusion_flag - πⱼ_prev) * neighbor_weights[i]
            end
        end
    end

    return BiodiversityObservationNetwork(pool[inclusion_flags], domain)
end
