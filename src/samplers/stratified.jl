"""
    SpatiallyStratified

`SpatiallyStratified` performs stratified random sampling over discrete
categories present in the domain. Each pool element belongs to a stratum given
by `domain[x]`. The number of draws allocated to each stratum is proportional
to the stratum size (via the multinomial distribution), and units are then sampled without replacement from each stratum.
"""
@kwdef struct SpatiallyStratified <: BONSampler
    num_nodes = _DEFAULT_NUM_NODES
end


"""
    _sample(sampler::SpatiallyStratified, domain; inclusion=nothing)

Draw a stratified random sample across unique values in `domain`.

Arguments:
- `sampler.num_nodes`: total number of units to sample across all strata
- `domain`: sampling domain; must support `getpool(domain)` and indexing `domain[x]`

Returns a `BiodiversityObservationNetwork`.
"""
function _sample(
    sampler::SpatiallyStratified,
    domain;
    inclusion = nothing
)
    
    pool = getpool(domain)
    vals = [domain[x] for x in pool]
    unique_vals = unique(vals)

    stratified_pool = [pool[findall(isequal(v), vals)] for v in unique_vals]

    area_weights = length.(stratified_pool) ./ length(pool)

    samples_per_stratum = rand(Multinomial(sampler.num_nodes, area_weights))

    nodes = vcat([SB.sample(stratified_pool[i], samples_per_stratum[i], replace=false) for i in eachindex(samples_per_stratum)]...)

    return nodes, domain[nodes]
end



