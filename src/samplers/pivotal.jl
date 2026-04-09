"""
    Pivotal <: BONSampler

Local Pivotal Method (LPM) for spatially balanced sampling.

Iteratively pairs nearby candidates and updates their inclusion probabilities
so that one tends toward selection and the other toward exclusion. Repeating
over the domain produces a spatially balanced sample that respects per-unit
inclusion probabilities.

# Fields
- `n::Int`: number of sites to select (default 50)

# References
- TODO Grafström, A. (2012).
"""
@kwdef struct Pivotal <: BONSampler
    n::Int = 50
end

supports_inclusion(::Pivotal) = true
guarantees_exact_n(::Pivotal) = true

"""
    _above_one_update!(inclusion, pool, i, j, complete_flags, selected_flag)

Apply the LPM update when the paired units have `πᵢ + πⱼ ≥ 1`.

One of the two units is set to 1 (selected) and the other is reduced to `πᵢ + πⱼ - 1`. The unit whose probability reaches 1 is marked both as included (`selected_flag`) and complete (`complete_flags`).
"""
function _above_one_update!(inclusion, i, j, complete_flags, selected_flag)
    πᵢ, πⱼ = inclusion[i], inclusion[j]
    add_idx, sub_idx = rand() < (1 - πⱼ)/(2 - πᵢ - πⱼ) ? (i, j) : (j, i)

    inclusion[add_idx] = 1     
    inclusion[sub_idx] = πᵢ + πⱼ - 1

    selected_flag[add_idx] = 1
    complete_flags[add_idx] = 1
end

"""
    _below_one_update!(inclusion, pool, i, j, complete_flags)

Apply the LPM update when the paired units have `πᵢ + πⱼ < 1`.

One of the two units is set to 0 (not selected) and the other is increased to
`πᵢ + πⱼ`. The unit whose probability reaches 0 is marked complete.
"""
function _below_one_update!(inclusion, i, j, complete_flags)
    πᵢ, πⱼ = inclusion[i], inclusion[j]
    add_idx, sub_idx = rand() < (πⱼ / (πᵢ + πⱼ)) ? (j, i) : (i, j)

    inclusion[add_idx] = πᵢ + πⱼ 
    inclusion[sub_idx] = 0
    complete_flags[sub_idx] = 1
end 

"""
    _apply_update_rule!(inclusion, pool, i, j, selected_flag, complete_flags)

Dispatch to the appropriate LPM update based on whether `πᵢ + πⱼ` is below or
at/above one.
"""
function _apply_update_rule!(inclusion, i, j, selected_flag, complete_flags)
    πᵢ, πⱼ = inclusion[i], inclusion[j]
    if πᵢ + πⱼ  < 1 
        _below_one_update!(inclusion, i, j, complete_flags)
    else
        _above_one_update!(inclusion, i, j, complete_flags, selected_flag)
    end 
end 


function _sample(rng::AbstractRNG, sampler::Pivotal, cpool::CandidatePool)
    N = cpool.n
    tree = _build_kdtree(cpool)

    # Inclusion probabilities scaled so they sum to n
    inclusion = cpool.inclusion .* sampler.n
    complete = zeros(Bool, N)
    selected = zeros(Bool, N)

    neighbor_map = _neighbor_map(tree, cpool.coordinates)
    
    while !all(complete)
        # Select a random incomplete index
        candidate_i = findall(!isone, complete)
        i = rand(rng, candidate_i)

        # find closest node to i that is not complete
        j_idx = findfirst(k->!complete[k], neighbor_map[i])
        
        if !isnothing(j_idx)
            j = neighbor_map[i][j_idx]
            _apply_update_rule!(inclusion, i, j, selected, complete)
        else 
            # if here, this is the last incomplete node 
            complete[i] = 1
        end     
    end

    return selected
end 

@testitem "We can use Pivotal with a BON" begin
    candidate_bon = sample(SimpleRandom(100), rand(30,20))

    bon = sample(Pivotal(50), candidate_bon)

    @test bon isa BiodiversityObservationNetwork
    @test first(bon) isa CartesianIndex
end

