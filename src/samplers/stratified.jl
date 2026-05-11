"""
    Stratified <: BONSampler

Stratified sampling. The first feature row of the `CandidatePool` is
treated as a discrete stratum label. Draws are allocated across strata via a
`Multinomial` draw on the stratum weights; within each stratum, sites are drawn
without replacement.

Recommended usage is to pass a single matrix/SDMLayer of discrete integer stratum labels.

# Fields
- `n::Int`: number of sites to select 
- `weights::Union{Vector, Missing}`: per-stratum weights. Index corresponds to label.
  `missing` (default) uses weights proportional to area.

# Notes
Does not guarantee exactly `n` selected sites when a Multinomial draw allocates
more draws to a stratum than it has candidates. In that case the stratum is
exhausted and the total falls short.

# Usage
```julia
weights = [0.1, 0.2, 0.5, 0.1, 0.1]
domain = rand(1:5, (30, 30))
sample(Stratified(weights = weights), domain)
```
"""
@kwdef struct Stratified <: BONSampler
    n::Int = 50
    weights::Union{Vector, Missing} = missing
end

requires_features(::Stratified) = true

function _sample(rng::AbstractRNG, sampler::Stratified, cpool::CandidatePool)
    labels = cpool.features[1, :]
    unique_labels = sort(unique(labels))

    strata = [findall(isequal(l), labels) for l in unique_labels]

    # Use area-based weights
    if ismissing(sampler.weights)
        w = length.(strata) ./ cpool.n
    else
        w = sampler.weights
        s = sum(w)
        s > 0 || throw(ArgumentError("Stratum weights must have a positive sum"))
        if s != 1 
            @info "Weights do not sum to one. Renormalize weights and proceeding"
            w ./= s
        end 
    end

    alloc = rand(rng, Distributions.Multinomial(sampler.n, w))

    selected = Int[]
    for (stratum_indices, k) in zip(strata, alloc)
        k == 0 && continue
        k = min(k, length(stratum_indices)) # in case allocated sites are greater than available number
        push!(selected, StatsBase.sample(rng, stratum_indices, k; replace=false)...)
    end

    return selected
end


@testitem "StratifiedRandom requires features" begin
    cp = CandidatePool(rand(10, 10))
    cp_nofeats = CandidatePool(cp.n, cp.keys, cp.coordinates, missing, cp.inclusion)
    @test_throws ArgumentError sample(Stratified(n=5), cp_nofeats)
end

@testitem "StratifiedRandom custom weights" begin
    cp = CandidatePool(rand(1:3, (20, 20)))
    result = sample(Stratified(n=20, weights=[0.4, 0.3, 0.3]), cp)
    @test result isa BiodiversityObservationNetwork
end
