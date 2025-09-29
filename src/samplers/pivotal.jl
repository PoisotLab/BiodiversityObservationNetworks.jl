"""
    Pivotal 

The Local Pivotal Method (2) from [Grafstrom2012SpaBal](@cite) is used for generating spatially balanced samples. 

`Pivotal` implements the Local Pivotal Method (LPM), which iteratively pairs nearby units and updates their inclusion probabilities so that, locally, one unit tends toward selection while the other tends toward non-selection. Repeating this over the domain produces a sample that is more spatially balanced than simple random sampling, while respecting per-unit inclusion probabilities when provided.

High-level algorithm:
- Select an unfinished unit `i`; then choose its nearest unfinished neighbor `j`.
- If `πᵢ + πⱼ < 1`, move probability mass so one unit goes to 0 and the other to `πᵢ+πⱼ`.
- Otherwise (if `πᵢ + πⱼ ≥ 1`), move probability mass so one unit is selected and the other retains the leftover probability `πᵢ + πⱼ - 1`.
- Repeat until all are complete.
"""
Base.@kwdef struct Pivotal{I<:Integer} <: BONSampler
    num_nodes::I = _DEFAULT_NUM_NODES
end 

"""
    _above_one_update!(inclusion, pool, i, j, complete_flags, inclusion_flags)

Apply the LPM update when the paired units have `πᵢ + πⱼ ≥ 1`.

One of the two units is set to 1 (selected) and the other is reduced to `πᵢ + πⱼ - 1`. The unit whose probability reaches 1 is marked both as included (`inclusion_flags`) and complete (`complete_flags`).
"""
function _above_one_update!(inclusion, pool, i, j, complete_flags, inclusion_flags)
    πᵢ, πⱼ = inclusion[pool[i]], inclusion[pool[j]] 
    add_idx, sub_idx = rand() < (1 - πⱼ)/(2 - πᵢ - πⱼ) ? (i, j) : (j, i)

    inclusion[pool[add_idx]] = 1     
    inclusion[pool[sub_idx]] = πᵢ + πⱼ - 1

    inclusion_flags[add_idx] = 1
    complete_flags[add_idx] = 1
end

"""
    _below_one_update!(inclusion, pool, i, j, complete_flags)

Apply the LPM update when the paired units have `πᵢ + πⱼ < 1`.

One of the two units is set to 0 (not selected) and the other is increased to
`πᵢ + πⱼ`. The unit whose probability reaches 0 is marked complete.
"""
function _below_one_update!(inclusion, pool, i, j, complete_flags)
    πᵢ, πⱼ = inclusion[pool[i]], inclusion[pool[j]] 
    add_idx, sub_idx = rand() < (πⱼ / (πᵢ + πⱼ)) ? (j, i) : (i, j)

    inclusion[pool[add_idx]] = πᵢ + πⱼ 
    inclusion[pool[sub_idx]] = 0

    complete_flags[sub_idx] = 1
end 

"""
    _apply_update_rule!(inclusion, pool, i, j, inclusion_flags, complete_flags)

Dispatch to the appropriate LPM update based on whether `πᵢ + πⱼ` is below or
at/above one.
"""
function _apply_update_rule!(inclusion, pool, i, j, inclusion_flags, complete_flags)
    πᵢ, πⱼ = inclusion[pool[i]], inclusion[pool[j]]
    if πᵢ + πⱼ  < 1 
        _below_one_update!(inclusion, pool, i, j, complete_flags)
    else
        _above_one_update!(inclusion, pool, i, j, complete_flags, inclusion_flags)
    end 
end 


"""
    _sample(sampler::Pivotal, domain; inclusion=nothing)

Draw a spatially balanced sample using the Local Pivotal Method.

Arguments:
- `sampler.num_nodes`: desired number of selected sites (also used to derive a
  uniform inclusion vector when `inclusion` is not provided)
- `domain`: sampling domain; must support `getpool(domain)` and `getnearestneighbors(domain)`
- `inclusion`: optional vector/array of inclusion probabilities indexed by pool items;
  if `nothing`, uniform probabilities are computed via `get_uniform_inclusion`.

Returns a `BiodiversityObservationNetwork` with nodes whose final inclusion indicators
are 1 after the LPM iterations.
"""
function _sample(sampler::Pivotal, domain; inclusion=nothing) 
    pool = getpool(domain)
    complete_flags = zeros(Bool, length(pool))

    inclusion = isnothing(inclusion) ? get_uniform_inclusion(sampler, domain) : inclusion
    
    inclusion_flags = zeros(Bool, length(pool))

    # Neighbor order for each unit, excluding itself
    # neighbor_order[i] => [index of unit closest to i, index of unit 2nd closest to i, ...]
    _, neighbor_order = getnearestneighbors(domain)

    while !all(complete_flags)
        candidate_i = findall(!isone, complete_flags)
        i = rand(candidate_i)

        # find closest node to i that is not complete
        j_idx = findfirst(k->!complete_flags[k], neighbor_order[i])
        if !isnothing(j_idx)
            j = neighbor_order[i][j_idx]
            _apply_update_rule!(inclusion, pool, i, j, inclusion_flags, complete_flags)
        else 
            complete_flags[i] = 1
        end 
    end
    return BiodiversityObservationNetwork(pool[findall(isone, inclusion_flags)], domain)
end 
