# Design of BiodiversityObservationNetworks.jl

This document provides an overview of how BiodiversityObservationNetworks.jl is designed, targeted at developers contributing to the package.

## Pipeline

For sampling sites, `sample(sampler, domain; mask, inclusion)` is the single public entry point. 

All domains are first transformed into a `CandidatePool(domain; ...)`, which stores information in a standardized format for sampling algorithms to act on.

Then, sampling is done in three stages. *First*, callign `validate(sampler, candidatepool)` to ensure the candidate pool meets the algorithm's requirements.

*Second*, the internal `_sample` method is called, which is where the logic defining each sampling algorithm is defined.  

*Finally*, the `_sample` method returns the indices of the selected sites within the `candidatepool.keys` vector, and this is used to construct the [`BiodiversityObservationNetwork`](@ref) type.

```
Domain -> CandidatePool -> Algorithm -> Selected Indices
                |                             |
                | ----------------------------|
                                              |         
                                              ↓
                                             BON
```

## Core Types

### [`CandidatePool`](@ref)

The [`CandidatePool`](@ref) type takes a domain representation and converts it to a standardized type for sampling algorithms to act on.  

| Field | Type | Notes |
|---|---|---|
| `n` | `Int` | number of candidates |
| `keys` | `Vector{K}` | original identifiers (e.g. `CartesianIndex{2}`) |
| `coordinates` | `Matrix` | `2 × n` spatial positions (either geospatial coordinates if using a SpeciesDistributionToolkit domain, or raster positions) |
| `features` | `Union{Matrix, Missing}` | `p × n` auxiliary variables associated with eac hsite |
| `inclusion` | `Vector` | inclusion probabiltiies, sums to 1 |

Inclusion always sums to 1. Samplers that need it to sum to `n` scale it internally as `π = cpool.inclusion .* sampler.n`.

### [`BONSampler`](@ref)

The [`BONSampler`](@ref) is an abstract supertype for all sampling algorithms. Every concrete sampler is subtype of [`BONSampler`](@ref) with at minimum `n::Int` as a field.

Several trait methods are defined to samplers, all of which default to `false`:

- `supports_inclusion(::BONSampler)`: can support non-uniform inclusion weights
- `supports_features(::BONSampler)`: can use auxiliary feature variables
- `requires_features(::BONSampler)`: throws error if features are missing
- `guarantees_exact_n(::BONSampler)`: will always return exactly `n` sites

### `BiodiversityObservationNetwork{K}`

The result of sampling. Carries selected sites, coordinates, features at those sites (or `Missing`), inclusion weights, and the sampler that produced it.

Can be passed back to `CandidatePool(bon)` for running multi-stage sampling with multiple algorithms.

## Implementing a New Sampler

Implementing a new sampling algorithm consists of three steps:

1. Defining a `@kwdef struct MySampler <: BONSampler` with `n::Int` and any additional parameters.
2. Implementing `_sample(rng::AbstractRNG, sampler::MySampler, cpool::CandidatePool) → Vector{Int}`, returning indices into `cpool.keys`.
3. Optionally override trait methods.

e.g.

```julia
@kwdef struct MySampler <: BONSampler
    n::Int = 50
end

# optional traits
guarantees_exact_n(::MySampler) = true

function _sample(rng, sampler::MySampler, cpool::CandidatePool)
    # return a Vector{Int} of length sampler.n, indexing into cpool.keys
end
```

The dispatch structure means these are the only methods required to implement a new sampler, and using `sample` as usual with any supported domain will work.

## Evaluation Metrics

All metrics are subtypes of `SamplingMetric`. The public function is `evaluate(metric, ...)`.

| Metric | Signature | Measures |
|---|---|---|
| `MoransI` | `(::MoransI, domain, bon)` | Spatial autocorrelation of the inclusion indicator. More negative = more spread |
| `VoronoiVariance` | `(::VoronoiVariance, domain, bon)` | Voronoi cell inclusion-weight variance; Closer to 0 = more spread |
| `JensenShannon` | `(::JensenShannon, layers, bon)` | Feature-space representativeness as JS divergence |


