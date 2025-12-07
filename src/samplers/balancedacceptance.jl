"""
    BalancedAcceptance <: BONSampler

Implements Balanced Acceptance Sampling (BAS) using Halton sequences.

# Fields
- `num_nodes::Int`: The number of sites to select.

# Description
BAS generates spatially balanced samples by mapping the domain to a Halton sequence.
If `inclusion` probabilities are provided, it uses a 3D Halton sequence where the 
third dimension acts as an acceptance threshold against the probability surface.

# References
- Robertson, B. L., et al. (2013). 
"""
@kwdef struct BalancedAcceptance <: BONSampler
    num_nodes = _DEFAULT_NUM_NODES
end 

_get_halton_value(base, offset, element) = haltonvalue(offset + element, base)
_halton(bases, seeds, step, dims) = [_get_halton_value(bases[i], seeds[i], step) for i in 1:dims] 


"""
    _rescale_node(domain, x::Real, y::Real)

Map unit-cube Halton coordinates `(x, y)` to integer raster indices in `domain`.
"""
function _rescale_node(domain, x::Real, y::Real) 
    x_scaled, y_scaled = Int.(ceil.(size(domain) .* [x, y]))
    return x_scaled, y_scaled
end    


"""
    _sample(sampler::BalancedAcceptance, domain; inclusion=nothing)

Generate a spatially balanced sample using BAS. With `inclusion`, perform 3D BAS
to respect per-cell probabilities; otherwise perform 2D BAS over the mask.
"""
function _sample(
    sampler::BalancedAcceptance,
    domain;
    inclusion = nothing
)   
    if !isnothing(inclusion)
        return _3d_bas(sampler, domain, inclusion)
    else
        return _2d_bas(sampler, domain)
    end
end


"""
    _3d_bas(sampler, domain, inclusion)

3D BAS using Halton bases `[2,3,5]`. A candidate `(i,j,z)` is accepted if the
cell is unmasked and `z < inclusion[i,j]`.
"""
function _3d_bas(sampler, domain, inclusion)
    seeds = rand(Int.(1e0:1e7), 3)
    bases = [2, 3, 5]
    attempt = 0
    nodes = []
    while length(nodes) < sampler.num_nodes
        i, j, z  = _halton(bases, seeds, attempt, 3)
        i, j = _rescale_node(domain, i,j)

        if !ismasked(domain, i,j) && z < inclusion[i,j] 
            push!(nodes, (i,j))
        end 
        attempt += 1
    end 
    return BiodiversityObservationNetwork(nodes, domain)
end 

"""
    _2d_bas(sampler, domain)

2D BAS using Halton bases `[2,3]` to generate spatially spread candidate cells,
accepting those that fall on unmasked locations until `num_nodes` are selected.
"""
function _2d_bas(sampler, domain)
    seeds = rand(Int.(1e0:1e7), 2)
    bases = [2, 3]
    attempt = 0
    nodes = []
    while length(nodes) < sampler.num_nodes
        i, j = _halton(bases, seeds, attempt, 2)
        i, j = _rescale_node(domain, i,j)
        
        if !ismasked(domain, i,j)
            push!(nodes, CartesianIndex(i,j))
        end 
        attempt += 1
    end
    return BiodiversityObservationNetwork(nodes, domain)
end 