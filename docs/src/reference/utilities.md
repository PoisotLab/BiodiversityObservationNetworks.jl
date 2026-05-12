# Utilities

Utilities for quantifying network quality and producing inputs (inclusion probabilities, rarity surfaces) useful when designing a [`BiodiversityObservationNetwork`](@ref).

## Spatial Balance

Two metrics evaluate how evenly selected sites cover the domain. Both are called via `evaluate(metric, domain, bon)`.

### [`MoransI`](@ref)

Spatial autocorrelation of the inclusion indicator (1 if sampled, 0 otherwise). Computes the normalized Moran's I statistic. Negative values indicate spatial inhibition (spread), which is desirable.

```@docs; canonical=false
MoransI
```

### [`VoronoiVariance`](@ref) 

Assigns every domain cell to its nearest selected site (Voronoi tessellation), then measures how far each cell's total inclusion weight deviates from 1. Defined as `(1/n) Σ(vᵢ - 1)²`; smaller is better.

```@docs; canonical=false
VoronoiVariance
```

## Environmental Representativeness

[`JensenShannon`](@ref) measures how well the selected sites represent the distribution of environmental predictors across the domain, using Jensen-Shannon divergence between per-feature empirical distributions. Returns average JS divergence across features, where smaller values mean the BON is more representative. 

```@docs; canonical=false
JensenShannon
```

## Climate Rarity

All rarity metrics are subtypes of `RarityMetric` and called via `evaluate(metric, layers)` or `evaluate(metric, layers, bon)` for BON-relative metrics. `layers` is a `Vector{<:Matrix}` (or `Vector{<:SDMLayer}` with SDT loaded). Output is a matrix of rarity scores.

| Metric | Requires BON | Description |
|---|---|---|
| `DistanceToMedian` | No | Euclidean distance in z-scored feature space to the per-feature median. Higher = more rare. |
| `MultivariateEnvironmentalSimilarity` | No | MESS surface. Per-cell minimum over features of an ECDF-derived similarity score. Negative indicates outside BON site range. |
| `DistanceToAnalogNode` | Yes | Distance in z-scored feature space to the nearest selected BON site. Higher = less represented. |
| `WithinRange` | Yes |  Whether each cell falls within the feature-space hyper-rectangle spanned by BON sites. |


```@docs; canonical=false
DistanceToMedian
MultivariateEnvironmentalSimilarity
DistanceToAnalogNode
WithinRange
```


