"""
    SimpleRandom
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