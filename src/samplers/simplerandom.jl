"""
    SimpleRandom <: BONSampler

Implements simple random sampling (SRS) and weighted random sampling.

# Fields
- `num_nodes::Int`: The number of sites to select.

# Description
Selects `num_nodes` locations from the domain uniformly at random without replacement. 
If an `inclusion` probability surface is provided, it performs weighted random sampling 
without replacement, where the weight of each candidate cell is proportional to its 
inclusion probability.

# Notes
While computationally efficient, Simple Random Sampling does not guarantee spatial 
balance and may result in clustering of sampling points.
"""
@kwdef struct SimpleRandom <: BONSampler
    num_nodes = _DEFAULT_NUM_NODES
end

function _sample(
    sampler::SimpleRandom,
    domain;
    inclusion = nothing
)
    pool = getpool(domain)

    nodes = isnothing(inclusion) ? SB.sample(pool, sampler.num_nodes, replace=false) : SB.wsample(pool, inclusion[pool], sampler.num_nodes, replace=false)

    return BiodiversityObservationNetwork(nodes, domain)
end 