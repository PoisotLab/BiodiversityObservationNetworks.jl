"""
    BalancedAcceptance <: BONSampler

Balanced Acceptance Sampling (BAS) using Halton sequences.

Generates spatially balanced samples by mapping Halton
sequences to the candidate coordinate space. When inclusion weights are
non-uniform, a third Halton dimension acts as a threshold for acceptance,
preferentially selecting higher-weighted candidates.

# Fields
- `n::Int`: number of sites to select (default 50)

# References
- Robertson, B. L., et al. (2013).
"""
@kwdef struct BalancedAcceptance <: BONSampler
    n::Int = 50
end

supports_inclusion(::BalancedAcceptance) = true
guarantees_exact_n(::BalancedAcceptance) = true

_get_halton_value(base, offset, element) = haltonvalue(offset + element, base)
_halton(bases, seeds, step, dims) = [_get_halton_value(bases[i], seeds[i], step) for i in 1:dims] 
_rescale_node(sz, x::Real, y::Real) = Int.(ceil.(sz .* [x, y]))

function _construct_explicit_mask_and_inclusion(cpool::CandidatePool, sz)
    mask, inclusion, key_index = zeros(Bool, sz), zeros(sz), zeros(Int, sz)

    for i in eachindex(cpool.keys)
        mask[cpool.keys[i]] = true
        inclusion[cpool.keys[i]] = cpool.inclusion[i]
        key_index[cpool.keys[i]] = i
    end
    return mask, inclusion, key_index
end

function _sample(rng::AbstractRNG, sampler::BalancedAcceptance, cpool::CandidatePool)
    unequal_inclusion = !all(≈(cpool.inclusion[1]), cpool.inclusion)
    if unequal_inclusion
        return _3d_bas(rng, sampler, cpool)
    else
        return _2d_bas(rng, sampler, cpool)
    end
end 

"""
    _2d_bas(sampler, domain)

2D BAS using Halton bases `[2,3]` to generate spatially spread candidate cells,
accepting those that fall on unmasked locations until `num_nodes` are selected.
"""
function _2d_bas(rng, sampler, cpool)
    seeds = rand(rng, Int.(1e0:1e7), 2)
    bases = [2, 3]
    attempt = 0
    selected = []
    xmax, ymax = maximum([x[1] for x in cpool.keys]), maximum([x[2] for x in cpool.keys])
    valid_sites, _, key_index = _construct_explicit_mask_and_inclusion(cpool, (xmax, ymax))

    while length(selected) < sampler.n
        i, j = _halton(bases, seeds, attempt, 2)
        i, j = _rescale_node((xmax, ymax), i, j)
        
        if valid_sites[i,j] 
            push!(selected, key_index[i,j])
        end 
        attempt += 1
    end
    return selected
end 

"""
    _3d_bas(rng, sampler, cpool)

3D BAS using Halton bases `[2,3,5]`. A candidate `(i,j,z)` is accepted if the
cell is unmasked and `z < inclusion[i,j]`.
"""
function _3d_bas(rng, sampler, cpool)
    seeds = rand(rng, Int.(1e0:1e7), 3)
    bases = [2, 3, 5]
    attempt = 0
    selected = []
    xmax, ymax = maximum([x[1] for x in cpool.keys]), maximum([x[2] for x in cpool.keys])
    
    valid_sites, inclusion, key_index = _construct_explicit_mask_and_inclusion(cpool, (xmax, ymax))

    while length(selected) < sampler.n
        i, j, z  = _halton(bases, seeds, attempt, 3)
        i, j = _rescale_node((xmax, ymax), i,j)
        
        if valid_sites[i,j]  && z < inclusion[i,j] 
            push!(selected, key_index[i,j])
        end 
        attempt += 1
    end 
    return selected
end 


@testitem "We can use Balanced Acceptance" begin
    bon = sample(BalancedAcceptance(), rand(30,20))
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can use Balanced Acceptance with custom inclusion" begin
    inclusion = rand(30,20)
    bon = sample(BalancedAcceptance(), rand(30,20), inclusion=inclusion)
    @test bon isa BiodiversityObservationNetwork
end

