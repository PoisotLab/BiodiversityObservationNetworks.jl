
# Public API {#Public-API}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.AdaptiveHotspot' href='#BiodiversityObservationNetworks.AdaptiveHotspot'><span class="jlbinding">BiodiversityObservationNetworks.AdaptiveHotspot</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
AdaptiveHotspot
```


Adaptive hotspot sampling prioritizes high-uncertainty regions while encouraging spatial diversity via a kernel-based criterion. This implementation follows ([Andrade-Pacheco _et al._, 2020](/bibliography#Andrade-Pacheco2020FinHot)): start at the maximum of the uncertainty surface and iteratively add locations that optimize a trade-off between local value and diversity with respect to already chosen sites using a Matérn covariance kernel.

Arguments:
- `num_nodes`: number of sites to select
  
- `scale` (`ρ`): range parameter of the Matérn kernel
  
- `smoothness` (`ν`): smoothness parameter of the Matérn kernel
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/adaptivehotspot.jl#L1-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.BONSampler' href='#BiodiversityObservationNetworks.BONSampler'><span class="jlbinding">BiodiversityObservationNetworks.BONSampler</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BONSampler
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/BiodiversityObservationNetworks.jl#L23-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.BalancedAcceptance' href='#BiodiversityObservationNetworks.BalancedAcceptance'><span class="jlbinding">BiodiversityObservationNetworks.BalancedAcceptance</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BalancedAcceptance
```


`BalancedAcceptance` implements Balanced Acceptance Sampling (BAS), which uses Halton sequences to provide spatially well-spread samples. When inclusion probabilities are supplied, a 3D Halton sequence is used: the first two dimensions spread locations spatially and the third acts as a random threshold for acceptance against the inclusion surface.

If inclusion probabilities are not provided, a 2D Halton sequence is used to generate a spatially balanced sample of size `num_nodes` while respecting raster masks.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/balancedacceptance.jl#L1-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.BiodiversityObservationNetwork' href='#BiodiversityObservationNetworks.BiodiversityObservationNetwork'><span class="jlbinding">BiodiversityObservationNetworks.BiodiversityObservationNetwork</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BiodiversityObservationNetwork
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/bon.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.CubeSampling' href='#BiodiversityObservationNetworks.CubeSampling'><span class="jlbinding">BiodiversityObservationNetworks.CubeSampling</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CubeSampling
```


`CubeSampling` implements the cube method of (???) for balanced sampling with respect to a set of auxiliary variables (features).

The algorithm proceeds in two phases:
- Flight phase: probabilities are iteratively moved towards 0/1 while preserving linear constraints on the mean of auxiliary variables in expectation (including fixed sample size).
  
- Landing phase: if fractional probabilities remain, an optimization step chooses a 0/1 sample that best matches the target constraints.
  

If `inclusion` is not provided, uniform inclusion probabilities are derived from `sampler.num_nodes` and the domain pool size. Returned nodes correspond to units with final probability 1.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L1-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.DistanceToMedian' href='#BiodiversityObservationNetworks.DistanceToMedian'><span class="jlbinding">BiodiversityObservationNetworks.DistanceToMedian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DistanceToMedian
```


Rarity score defined as Euclidean distance in feature space to the per-feature median across the raster stack. Optionally, features can be PCA-transformed prior to z-scoring.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/rarity.jl#L27-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.GeneralizedRandomTesselated' href='#BiodiversityObservationNetworks.GeneralizedRandomTesselated'><span class="jlbinding">BiodiversityObservationNetworks.GeneralizedRandomTesselated</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GeneralizedRandomTesselated
```


`GeneralizedRandomTesselated` is a type of [`BONSampler`](/reference/api#BiodiversityObservationNetworks.BONSampler) for generating [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork)s with spatial spreading.

GRTS was initially proposed in ([Stevens and Olsen, 2004](/bibliography#Stevens2004SpaBal)).

_Arguments_:
- `num_nodes`: the number of sites to select
  
- `grid_size`: if being used on a polygon, the dimensions of the grid used to cover the extent. GRTS sampling uses discrete Cartesian indices
  

GRTS represents each cell of a rasterized version of the sampling domain using an address, where the address of each cell is represented as a `D`-digit base-4 number. 

The value of `D` depends on the size of the raster. GRTS works by recursively splitting the rasterized domain into quadrants, and those quadrants into further quadrants, until a single pixel is reached. At each step, each quadrant is randomly labeled with a number 1, 2, 3, or 4 (without replacement, so all values are used). 

The addresses are then sorted numerically, and the `num_nodes` smallest values are chosen.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/grts.jl#L1-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.MoransI' href='#BiodiversityObservationNetworks.MoransI'><span class="jlbinding">BiodiversityObservationNetworks.MoransI</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MoransI
```


`MoransI` is a measure of spatial balance proposed by ([Tillé _et al._, 2018](/bibliography#Tille2018MeaSpa)).

Conceptually, the idea is to use [Moran&#39;s I](https://en.wikipedia.org/wiki/Moran%27s_I), a measure of spatial autocorrelation, on an indicator variable $\delta_i$, which is $1$ if unit $i$ is included in the sample, and 0 otherwise. 

If $\delta$ has a negative value, this means included samples have negative autocorrelation and therefore are more spread out than if generated at random. 

In principle, Moran&#39;s I should exist on the interval $[-1,1]$, but the original definition uses a renormalization that doesn&#39;t always guarantee this, so ([Tillé _et al._, 2018](/bibliography#Tille2018MeaSpa)) introduce a slight variation, called $I_C$ so this always holds.

Note that the performance of this algorithm scales as $O(N^3)$, where $N$ is the number of all possible candidate locations, and therefore for large rasters this method is unlikely to be performant. A faster approximation could be made by taking a smaller random sample of the candidate points. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/spatialbalance.jl#L15-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.MultivariateEnvironmentalSimilarity' href='#BiodiversityObservationNetworks.MultivariateEnvironmentalSimilarity'><span class="jlbinding">BiodiversityObservationNetworks.MultivariateEnvironmentalSimilarity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MultivariateEnvironmentalSimilarity
```


Multivariate Environmental Similarity Surface (MESS). For each cell, compute the minimum over features of a per-feature similarity score derived from the ECDF of the training distribution, following the standard MESS definition.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/rarity.jl#L57-L63" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.Pivotal' href='#BiodiversityObservationNetworks.Pivotal'><span class="jlbinding">BiodiversityObservationNetworks.Pivotal</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Pivotal
```


The Local Pivotal Method (2) from ([Grafström _et al._, 2012](/bibliography#Grafstrom2012SpaBal)) is used for generating spatially balanced samples. 

`Pivotal` implements the Local Pivotal Method (LPM), which iteratively pairs nearby units and updates their inclusion probabilities so that, locally, one unit tends toward selection while the other tends toward non-selection. Repeating this over the domain produces a sample that is more spatially balanced than simple random sampling, while respecting per-unit inclusion probabilities when provided.

High-level algorithm:
- Select an unfinished unit `i`; then choose its nearest unfinished neighbor `j`.
  
- If `πᵢ + πⱼ < 1`, move probability mass so one unit goes to 0 and the other to `πᵢ+πⱼ`.
  
- Otherwise (if `πᵢ + πⱼ ≥ 1`), move probability mass so one unit is selected and the other retains the leftover probability `πᵢ + πⱼ - 1`.
  
- Repeat until all are complete.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/pivotal.jl#L1-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.SimpleRandom' href='#BiodiversityObservationNetworks.SimpleRandom'><span class="jlbinding">BiodiversityObservationNetworks.SimpleRandom</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SimpleRandom
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/simplerandom.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.SpatiallyCorrelatedPoisson' href='#BiodiversityObservationNetworks.SpatiallyCorrelatedPoisson'><span class="jlbinding">BiodiversityObservationNetworks.SpatiallyCorrelatedPoisson</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SpatiallyCorrelatedPoisson
```


Implements the Spatially Correlated Poisson Sampling algorithm from Grafström (2012), Journal of Statistical Planning and Inference. doi:10.1016/j.jspi.2011.07.003

Iterate through all potential nodes i. Select it w.p. πᵢ, and adjust all πⱼ for all subsequent nodes j to maintain expected sample size.

I think the issue w/ samples occasionally being off-by-one is the last sample is rejected and the remaining inclusion probability has nowhere left to go. Same for Pivotal.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/spatiallycorrelatedpoisson.jl#L1-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.SpatiallyStratified' href='#BiodiversityObservationNetworks.SpatiallyStratified'><span class="jlbinding">BiodiversityObservationNetworks.SpatiallyStratified</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SpatiallyStratified
```


`SpatiallyStratified` performs stratified random sampling over discrete categories present in the domain. Each pool element belongs to a stratum given by `domain[x]`. The number of draws allocated to each stratum is proportional to the stratum size (via the multinomial distribution), and units are then sampled without replacement from each stratum.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/stratified.jl#L1-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.VoronoiVariance' href='#BiodiversityObservationNetworks.VoronoiVariance'><span class="jlbinding">BiodiversityObservationNetworks.VoronoiVariance</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
VoronoiVariance
```


The `VoronoiVariance` method for characterizing the spatial balance of a sample is based on the initial method proposed by ([Stevens and Olsen, 2004](/bibliography#Stevens2004SpaBal)), and then extended by ([Grafström, 2012](/bibliography#Grafstrom2012SpaCor)).

For a given [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork) `bon`, the [Voronoi tesselation](https://en.wikipedia.org/wiki/Voronoi_diagram) splits the plane into a series of polygons, where the `i`-th polygon consists of all points in the plane whose nearest node in `bon` is the `i`-th node.

These polygons can then be used to assess the spatial balance of a sample. 

In an _ideally_ balanced sample, the sum of the inclusion probabilities across each polygon `i` would equal 1, because in expectation exactly one unit would be sampled in that region. 

If we define $v_i$ as the total inclusion probability across all elements of the population in Voronoi polygon `i`, i.e.

$$v_i = \sum_{j \in i} \pi_j $$

then we can assess the spatial balance of a sample by measuring the distance of $v_i$ from 1 for each polygon.  ([Grafström, 2012](/bibliography#Grafstrom2012SpaCor)) proposes the metric `B`, defined as 

$$B = \frac{1}{n} \sum_{i=1}^n (v_i - 1)^2$$

to measure spatial balance, where _smaller values indicate better spatial balance_.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/spatialbalance.jl#L94-L128" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.extent-Tuple{PolygonDomain{<:SimpleSDMPolygons.AbstractGeometry}}' href='#BiodiversityObservationNetworks.extent-Tuple{PolygonDomain{<:SimpleSDMPolygons.AbstractGeometry}}'><span class="jlbinding">BiodiversityObservationNetworks.extent</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
extent(p::PolygonDomain{T<:SimpleSDMPolygons.AbstractGeometry})
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/polygon.jl#L5-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.extent-Tuple{RasterDomain{<:Matrix}}' href='#BiodiversityObservationNetworks.extent-Tuple{RasterDomain{<:Matrix}}'><span class="jlbinding">BiodiversityObservationNetworks.extent</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
extent
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L49-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.jensenshannon-Tuple{RasterStack, BiodiversityObservationNetwork}' href='#BiodiversityObservationNetworks.jensenshannon-Tuple{RasterStack, BiodiversityObservationNetwork}'><span class="jlbinding">BiodiversityObservationNetworks.jensenshannon</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
jensenshannon
```


The [Jensen-Shannon Divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence) is a method for measuring the distance between two probability distributions.

This method provides a comparison between the distribution of environmental variables in set of `SDMLayers`, `layers`, to the values of those variables at the sites within a [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork) `bon`. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/distances.jl#L28-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.rarity-Tuple{DistanceToAnalogNode, BiodiversityObservationNetwork, RasterStack}' href='#BiodiversityObservationNetworks.rarity-Tuple{DistanceToAnalogNode, BiodiversityObservationNetwork, RasterStack}'><span class="jlbinding">BiodiversityObservationNetworks.rarity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rarity(::DistanceToAnalogNode, bon, layers; pca=false)
```


For each cell, compute the distance in z-scored feature space to the nearest selected BON node in `layers`. Optionally apply a shared PCA transform first.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/rarity.jl#L100-L105" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.rarity-Tuple{WithinRange, BiodiversityObservationNetwork, RasterStack}' href='#BiodiversityObservationNetworks.rarity-Tuple{WithinRange, BiodiversityObservationNetwork, RasterStack}'><span class="jlbinding">BiodiversityObservationNetworks.rarity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rarity(::WithinRange, bon, layers)
```


Boolean rarity surface indicating whether each cell lies within the hyper- rectangle spanned by the per-feature minima and maxima of the BON nodes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/rarity.jl#L150-L155" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.sample-Tuple{}' href='#BiodiversityObservationNetworks.sample-Tuple{}'><span class="jlbinding">BiodiversityObservationNetworks.sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sample
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/sample.jl#L14-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.spatialbalance-Tuple{MoransI, BiodiversityObservationNetwork, RasterDomain}' href='#BiodiversityObservationNetworks.spatialbalance-Tuple{MoransI, BiodiversityObservationNetwork, RasterDomain}'><span class="jlbinding">BiodiversityObservationNetworks.spatialbalance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spatialbalance(::MoransI, raster::SDMLayer, bon::BiodiversityObservationNetwork)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/spatialbalance.jl#L58-L60" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.spatialbalance-Tuple{VoronoiVariance, BiodiversityObservationNetwork, PolygonDomain}' href='#BiodiversityObservationNetworks.spatialbalance-Tuple{VoronoiVariance, BiodiversityObservationNetwork, PolygonDomain}'><span class="jlbinding">BiodiversityObservationNetworks.spatialbalance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spatialbalance(::Type{VoronoiVariance}, bon::BiodiversityObservationNetwork, geom)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/spatialbalance.jl#L133-L135" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.spatialbalance-Union{Tuple{T}, Tuple{Type{T}, Any, Any}} where T' href='#BiodiversityObservationNetworks.spatialbalance-Union{Tuple{T}, Tuple{Type{T}, Any, Any}} where T'><span class="jlbinding">BiodiversityObservationNetworks.spatialbalance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spatialbalance
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/spatialbalance.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.voronoi-Tuple{Any, Any}' href='#BiodiversityObservationNetworks.voronoi-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks.voronoi</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
voronoi(bon, domain)
```


Construct Voronoi polygons for the nodes in a `BiodiversityObservationNetwork` within the given domain. The domain is coerced to a `PolygonDomain` when needed. Output polygons are clipped to the domain extent.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/voronoi.jl#L1-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>


# Private API {#Private-API}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.RarityMetric' href='#BiodiversityObservationNetworks.RarityMetric'><span class="jlbinding">BiodiversityObservationNetworks.RarityMetric</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RarityMetric
```


Abstract type encompassing all methods for computing environmental rarity.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/rarity.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._2d_bas-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._2d_bas-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._2d_bas</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_2d_bas(sampler, domain)
```


2D BAS using Halton bases `[2,3]` to generate spatially spread candidate cells, accepting those that fall on unmasked locations until `num_nodes` are selected.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/balancedacceptance.jl#L72-L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._3d_bas-Tuple{Any, Any, Any}' href='#BiodiversityObservationNetworks._3d_bas-Tuple{Any, Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._3d_bas</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_3d_bas(sampler, domain, inclusion)
```


3D BAS using Halton bases `[2,3,5]`. A candidate `(i,j,z)` is accepted if the cell is unmasked and `z < inclusion[i,j]`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/balancedacceptance.jl#L49-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._above_one_update!-NTuple{6, Any}' href='#BiodiversityObservationNetworks._above_one_update!-NTuple{6, Any}'><span class="jlbinding">BiodiversityObservationNetworks._above_one_update!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_above_one_update!(inclusion, pool, i, j, complete_flags, inclusion_flags)
```


Apply the LPM update when the paired units have `πᵢ + πⱼ ≥ 1`.

One of the two units is set to 1 (selected) and the other is reduced to `πᵢ + πⱼ - 1`. The unit whose probability reaches 1 is marked both as included (`inclusion_flags`) and complete (`complete_flags`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/pivotal.jl#L18-L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._apply_update_rule!-NTuple{6, Any}' href='#BiodiversityObservationNetworks._apply_update_rule!-NTuple{6, Any}'><span class="jlbinding">BiodiversityObservationNetworks._apply_update_rule!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_apply_update_rule!(inclusion, pool, i, j, inclusion_flags, complete_flags)
```


Dispatch to the appropriate LPM update based on whether `πᵢ + πⱼ` is below or at/above one.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/pivotal.jl#L54-L59" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._below_one_update!-NTuple{5, Any}' href='#BiodiversityObservationNetworks._below_one_update!-NTuple{5, Any}'><span class="jlbinding">BiodiversityObservationNetworks._below_one_update!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_below_one_update!(inclusion, pool, i, j, complete_flags)
```


Apply the LPM update when the paired units have `πᵢ + πⱼ < 1`.

One of the two units is set to 0 (not selected) and the other is increased to `πᵢ + πⱼ`. The unit whose probability reaches 0 is marked complete.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/pivotal.jl#L36-L43" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._cube_flight_phase-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._cube_flight_phase-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._cube_flight_phase</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_cube_flight_phase(πₖ, x)
```


Run the flight phase of the cube method.

`πₖ` are current inclusion probabilities; `x` (auxiliary matrix) is augmented with a first row of `πₖ'` so sample size fixed. The method repeatedly finds a direction in the null space of the constraint matrix and pushes a small subset of probabilities to 0 or 1 while preserving balances in expectation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L89-L97" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._cube_landing_phase-Tuple{Any, Any, Any}' href='#BiodiversityObservationNetworks._cube_landing_phase-Tuple{Any, Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._cube_landing_phase</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_cube_landing_phase(pikstar, πₖ, x)
```


Run the landing phase if some probabilities are still fractional after flight. Formulate a small linear program over the fractional units to select a 0/1 sample whose auxiliary totals match `pikstar` in expectation with minimal cost `C(s)` (per Deville &amp; Tillé, 2004).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L191-L198" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._get_address_length-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._get_address_length-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._get_address_length</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_get_address_length
```


Returns the number of digits in a [`GeneralizedRandomTessellatedStratified`](@ref) address for a specific geometry, computed based on the raster dimensions. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/grts.jl#L80-L85" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._get_addresses-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._get_addresses-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._get_addresses</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_get_addresses()
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/grts.jl#L88-L90" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._h-Tuple{Any}' href='#BiodiversityObservationNetworks._h-Tuple{Any}'><span class="jlbinding">BiodiversityObservationNetworks._h</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_h(K)
```


Entropy-like diversity score based on the log-determinant of the kernel matrix `K`. Larger values encourage selecting points that are diverse with respect to previously chosen sites.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/adaptivehotspot.jl#L89-L95" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._matérn-Tuple{Any, Any, Any}' href='#BiodiversityObservationNetworks._matérn-Tuple{Any, Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._matérn</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_matérn(d, ρ, ν)
```


Matérn covariance kernel evaluated at distance `d`, range `ρ`, smoothness `ν`. Normalized so that `_matérn(0, ρ, ν) == 1`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/adaptivehotspot.jl#L76-L81" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._pick_nodes-Tuple{Any, Any, Any}' href='#BiodiversityObservationNetworks._pick_nodes-Tuple{Any, Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._pick_nodes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_pick_nodes(sampler, raster, addresses)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/grts.jl#L104-L106" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._quadrant_fill!-Tuple{Any}' href='#BiodiversityObservationNetworks._quadrant_fill!-Tuple{Any}'><span class="jlbinding">BiodiversityObservationNetworks._quadrant_fill!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_quadrant_fill!(mat)
```


Takes a matrix `mat` and splits it into quadrants randomly labeled one through four. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/grts.jl#L44-L48" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._quadrant_split!-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._quadrant_split!-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._quadrant_split!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_quadrant_split!(mat, grid_size)
```


Splits a matrix `mat` into nested quadrants, where the side-length of a submatrix to be split into quadrants is given by `grid_size`.  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/grts.jl#L61-L66" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._redistribute_masked_weight!-Tuple{Any, RasterDomain{<:Matrix}}' href='#BiodiversityObservationNetworks._redistribute_masked_weight!-Tuple{Any, RasterDomain{<:Matrix}}'><span class="jlbinding">BiodiversityObservationNetworks._redistribute_masked_weight!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_redistribute_masked_weight!(domain, inclusion::RasterDomain{<:Matrix})
```


Matrix-backed variant of masked weight redistribution.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/inclusion.jl#L48-L52" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._redistribute_masked_weight!-Tuple{Any, RasterDomain{<:SimpleSDMLayers.SDMLayer}}' href='#BiodiversityObservationNetworks._redistribute_masked_weight!-Tuple{Any, RasterDomain{<:SimpleSDMLayers.SDMLayer}}'><span class="jlbinding">BiodiversityObservationNetworks._redistribute_masked_weight!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_redistribute_masked_weight!(domain, inclusion::RasterDomain{<:SDMLayer})
```


Redistribute any probability mass assigned to masked-out cells uniformly across valid indices in `domain`. Mutates `inclusion` in-place.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/inclusion.jl#L32-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._rescale_node-Tuple{Any, Real, Real}' href='#BiodiversityObservationNetworks._rescale_node-Tuple{Any, Real, Real}'><span class="jlbinding">BiodiversityObservationNetworks._rescale_node</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_rescale_node(domain, x::Real, y::Real)
```


Map unit-cube Halton coordinates `(x, y)` to integer raster indices in `domain`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/balancedacceptance.jl#L19-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._sample-Tuple{BalancedAcceptance, Any}' href='#BiodiversityObservationNetworks._sample-Tuple{BalancedAcceptance, Any}'><span class="jlbinding">BiodiversityObservationNetworks._sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_sample(sampler::BalancedAcceptance, domain; inclusion=nothing)
```


Generate a spatially balanced sample using BAS. With `inclusion`, perform 3D BAS to respect per-cell probabilities; otherwise perform 2D BAS over the mask.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/balancedacceptance.jl#L30-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._sample-Tuple{CubeSampling, Any}' href='#BiodiversityObservationNetworks._sample-Tuple{CubeSampling, Any}'><span class="jlbinding">BiodiversityObservationNetworks._sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_sample(sampler::CubeSampling, domain; inclusion=nothing)
```


Draw a sample using the cube method, balancing means of auxiliary variables (`getfeatures(domain)`) while achieving the desired sample size in expectation.

Arguments:
- `sampler.num_nodes`: desired number of selected sites
  
- `domain`: sampling domain; must support `getpool(domain)` and `getfeatures(domain)`
  
- `inclusion`: optional vector/array of inclusion probabilities indexed by pool items
  

Returns a `BiodiversityObservationNetwork`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L21-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._sample-Tuple{Pivotal, Any}' href='#BiodiversityObservationNetworks._sample-Tuple{Pivotal, Any}'><span class="jlbinding">BiodiversityObservationNetworks._sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_sample(sampler::Pivotal, domain; inclusion=nothing)
```


Draw a spatially balanced sample using the Local Pivotal Method.

Arguments:
- `sampler.num_nodes`: desired number of selected sites (also used to derive a uniform inclusion vector when `inclusion` is not provided)
  
- `domain`: sampling domain; must support `getpool(domain)` and `getnearestneighbors(domain)`
  
- `inclusion`: optional vector/array of inclusion probabilities indexed by pool items; if `nothing`, uniform probabilities are computed via `get_uniform_inclusion`.
  

Returns a `BiodiversityObservationNetwork` with nodes whose final inclusion indicators are 1 after the LPM iterations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/pivotal.jl#L70-L84" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._sample-Tuple{SpatiallyCorrelatedPoisson, Any}' href='#BiodiversityObservationNetworks._sample-Tuple{SpatiallyCorrelatedPoisson, Any}'><span class="jlbinding">BiodiversityObservationNetworks._sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_sample(sampler::SpatiallyCorrelatedPoisson, domain, inclusion_probs)
```


Perform Spatially Correlated Poisson sampling following Grafström (2012).

**Arguments**
- `sampler::SpatiallyCorrelatedPoisson`: The sampling strategy object.
  
- `domain`: The spatial domain (supports `size(domain)` and distance calculation).
  
- `inclusion_probs`: Initial inclusion probabilities for each unit.       Defaults to uniform: `(number_of_nodes / N)` for all units.
  

**Returns**
- A [`BiodiversityObservationNetwork.jl`](@ref) containing units included in the sample.
  

**Example Usage**

`sample(SpatiallyCorrelatedPoisson(), zeros(30,20))`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/spatiallycorrelatedpoisson.jl#L22-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._sample-Tuple{SpatiallyStratified, Any}' href='#BiodiversityObservationNetworks._sample-Tuple{SpatiallyStratified, Any}'><span class="jlbinding">BiodiversityObservationNetworks._sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_sample(sampler::SpatiallyStratified, domain; inclusion=nothing)
```


Draw a stratified random sample across unique values in `domain`.

Arguments:
- `sampler.num_nodes`: total number of units to sample across all strata
  
- `domain`: sampling domain; must support `getpool(domain)` and indexing `domain[x]`
  

Returns a `BiodiversityObservationNetwork`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/stratified.jl#L14-L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._sort_features_by_mahalanobis-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._sort_features_by_mahalanobis-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._sort_features_by_mahalanobis</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_sort_features_by_mahalanobis(features, inclusion)
```


Order units to stabilize the flight phase by spreading early decisions across feature space. Units are sorted by Mahalanobis distance from the inclusion- weighted mean feature vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L68-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._standardize-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._standardize-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._standardize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_standardize
```


Standardizes the values of a matrix of predictors across the entire population `Xfull`, and a set of predictors associated with the sampled sites, `Xsampled` by scaling each predictor to [0,1].

`Xsampled` is standardized based on the minimum and maximum values of each predictor across the population, so both matrices a return on the same scale.

_Arguments_:
- `Xfull`: an `n` x `d` matrix, where `n` is the size of the population, and `d` is the number of predictors
  
- `Xsampled`: an `m` x `d` matrix, where `m` &lt; `n` is the size of the sample
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/distances.jl#L1-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks._update_psi-Tuple{Any, Any}' href='#BiodiversityObservationNetworks._update_psi-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks._update_psi</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_update_psi(Ψ, B)
```


Given current working vector `Ψ` and constraint block `B`, compute a feasible update along a null-space direction `u` by maximizing step sizes `λ₁, λ₂` that keep probabilities within [0,1].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L147-L153" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_inclusion-Tuple{Any, Any, Nothing}' href='#BiodiversityObservationNetworks.convert_inclusion-Tuple{Any, Any, Nothing}'><span class="jlbinding">BiodiversityObservationNetworks.convert_inclusion</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_inclusion(sampler, domain, inclusion; kwargs...)
```


Coerce an `inclusion` specification to a representation compatible with `domain`. Also ensures totals sum to `sampler.num_nodes` (renormalizing if needed) and redistributes probability assigned to masked-out cells back to valid cells. Returns `nothing` when `inclusion` is `nothing`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/inclusion.jl#L89-L96" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, AbstractMatrix}' href='#BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, AbstractMatrix}'><span class="jlbinding">BiodiversityObservationNetworks.convert_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_mask(domain::RasterDomain, mask::AbstractMatrix; kwargs...)
```


Convert a boolean matrix mask to a `RasterDomain`-aligned mask.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L20-L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, PolygonDomain}' href='#BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, PolygonDomain}'><span class="jlbinding">BiodiversityObservationNetworks.convert_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_mask(domain::RasterDomain, mask::PolygonDomain; kwargs...)
```


Rasterize a polygon mask into the grid of `domain` and return a `RasterDomain`-aligned mask.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L71-L76" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, SimpleSDMLayers.SDMLayer}' href='#BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, SimpleSDMLayers.SDMLayer}'><span class="jlbinding">BiodiversityObservationNetworks.convert_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_mask(domain::RasterDomain, mask::SDMLayer; kwargs...)
```


Convert an `SDMLayer` to a `RasterDomain` mask with aligned indices.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L30-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, SimpleSDMPolygons.AbstractGeometry}' href='#BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain, SimpleSDMPolygons.AbstractGeometry}'><span class="jlbinding">BiodiversityObservationNetworks.convert_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_mask(domain::RasterDomain, mask::SimpleSDMPolygons.AbstractGeometry)
```


Coerce a polygon geometry to a `RasterDomain` mask by rasterizing the polygon.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L62-L66" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain{<:SimpleSDMLayers.SDMLayer}, RasterDomain{<:SimpleSDMLayers.SDMLayer}}' href='#BiodiversityObservationNetworks.convert_mask-Tuple{RasterDomain{<:SimpleSDMLayers.SDMLayer}, RasterDomain{<:SimpleSDMLayers.SDMLayer}}'><span class="jlbinding">BiodiversityObservationNetworks.convert_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_mask(domain::RasterDomain{<:SDMLayer}, mask::RasterDomain{<:SDMLayer}; kwargs...)
```


Validate size, extent, and CRS when both domain and mask are SDMLayer-backed.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L39-L43" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_mask-Union{Tuple{T}, Tuple{RasterDomain{T}, RasterDomain{<:Matrix}}} where T<:Union{var"#s23", var"#s24"} where {var"#s23"<:SimpleSDMLayers.SDMLayer, var"#s24"<:Matrix}' href='#BiodiversityObservationNetworks.convert_mask-Union{Tuple{T}, Tuple{RasterDomain{T}, RasterDomain{<:Matrix}}} where T<:Union{var"#s23", var"#s24"} where {var"#s23"<:SimpleSDMLayers.SDMLayer, var"#s24"<:Matrix}'><span class="jlbinding">BiodiversityObservationNetworks.convert_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_mask(domain::RasterDomain{T}, mask::RasterDomain{<:Matrix}; kwargs...) where T
```


Validate size when the mask is matrix-backed.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L52-L56" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_node-Tuple{RasterDomain{<:Matrix}, Any}' href='#BiodiversityObservationNetworks.convert_node-Tuple{RasterDomain{<:Matrix}, Any}'><span class="jlbinding">BiodiversityObservationNetworks.convert_node</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_node
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L63-L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.convert_nodes-Tuple{RasterDomain, Any}' href='#BiodiversityObservationNetworks.convert_nodes-Tuple{RasterDomain, Any}'><span class="jlbinding">BiodiversityObservationNetworks.convert_nodes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_nodes
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L69-L71" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.crs-Tuple{RasterDomain{<:Matrix}}' href='#BiodiversityObservationNetworks.crs-Tuple{RasterDomain{<:Matrix}}'><span class="jlbinding">BiodiversityObservationNetworks.crs</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
crs
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L56-L58" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.get_uniform_inclusion-Tuple{Any, BiodiversityObservationNetwork}' href='#BiodiversityObservationNetworks.get_uniform_inclusion-Tuple{Any, BiodiversityObservationNetwork}'><span class="jlbinding">BiodiversityObservationNetworks.get_uniform_inclusion</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_uniform_inclusion(sampler, bon::BiodiversityObservationNetwork)
```


Construct a uniform inclusion vector for a BON-like domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/inclusion.jl#L18-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.get_uniform_inclusion-Union{Tuple{T}, Tuple{Any, T}} where T<:Union{RasterDomain, RasterStack}' href='#BiodiversityObservationNetworks.get_uniform_inclusion-Union{Tuple{T}, Tuple{Any, T}} where T<:Union{RasterDomain, RasterStack}'><span class="jlbinding">BiodiversityObservationNetworks.get_uniform_inclusion</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_uniform_inclusion(sampler, domain)
```


Construct a uniform inclusion surface for `domain` such that the sum of inclusion probabilities equals `sampler.num_nodes`.

Returns a `RasterDomain` when given a raster-like domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/inclusion.jl#L1-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.getcoordinates-Tuple{BiodiversityObservationNetwork}' href='#BiodiversityObservationNetworks.getcoordinates-Tuple{BiodiversityObservationNetwork}'><span class="jlbinding">BiodiversityObservationNetworks.getcoordinates</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getcoordinates
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/bon.jl#L37-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.getcoordinates-Tuple{RasterDomain{<:Matrix}}' href='#BiodiversityObservationNetworks.getcoordinates-Tuple{RasterDomain{<:Matrix}}'><span class="jlbinding">BiodiversityObservationNetworks.getcoordinates</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getcoordinates
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L37-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.getfeatures-Tuple{BiodiversityObservationNetwork}' href='#BiodiversityObservationNetworks.getfeatures-Tuple{BiodiversityObservationNetwork}'><span class="jlbinding">BiodiversityObservationNetworks.getfeatures</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getfeatures
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/bon.jl#L32-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.getfeatures-Tuple{RasterDomain}' href='#BiodiversityObservationNetworks.getfeatures-Tuple{RasterDomain}'><span class="jlbinding">BiodiversityObservationNetworks.getfeatures</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getfeatures
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L29-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.getpool-Tuple{BiodiversityObservationNetwork}' href='#BiodiversityObservationNetworks.getpool-Tuple{BiodiversityObservationNetwork}'><span class="jlbinding">BiodiversityObservationNetworks.getpool</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getpool
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/bon.jl#L27-L29" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.getpool-Tuple{RasterDomain}' href='#BiodiversityObservationNetworks.getpool-Tuple{RasterDomain}'><span class="jlbinding">BiodiversityObservationNetworks.getpool</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getpool
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L23-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.ismasked-Tuple{RasterDomain{<:SimpleSDMLayers.SDMLayer}, Real, Real}' href='#BiodiversityObservationNetworks.ismasked-Tuple{RasterDomain{<:SimpleSDMLayers.SDMLayer}, Real, Real}'><span class="jlbinding">BiodiversityObservationNetworks.ismasked</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ismasked
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L90-L92" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.mask!-Tuple{Any, Nothing}' href='#BiodiversityObservationNetworks.mask!-Tuple{Any, Nothing}'><span class="jlbinding">BiodiversityObservationNetworks.mask!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mask!(domain, mask)
```


Apply `mask` to `domain` in-place. Overloads handle `RasterDomain`s backed by `SDMLayer` or `Matrix`, and propagate to all rasters in a `RasterStack`. Polygon masks are supported for SDMLayer-backed rasters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/mask.jl#L88-L94" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.normalize_inclusion!-Tuple{Any, RasterDomain{<:SimpleSDMLayers.SDMLayer}}' href='#BiodiversityObservationNetworks.normalize_inclusion!-Tuple{Any, RasterDomain{<:SimpleSDMLayers.SDMLayer}}'><span class="jlbinding">BiodiversityObservationNetworks.normalize_inclusion!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
normalize_inclusion!(sampler, inclusion)
```


Scale `inclusion` so that its total mass equals `sampler.num_nodes`. Overloads exist for SDMLayer-backed rasters, matrix-backed rasters, and vectors. Mutates `inclusion` in-place.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/inclusion.jl#L66-L72" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.rescale_node-Tuple{RasterDomain{<:Matrix}, Real, Real}' href='#BiodiversityObservationNetworks.rescale_node-Tuple{RasterDomain{<:Matrix}, Real, Real}'><span class="jlbinding">BiodiversityObservationNetworks.rescale_node</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rescale_node
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/domains/raster.jl#L75-L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.tilt-Tuple{Any, Any}' href='#BiodiversityObservationNetworks.tilt-Tuple{Any, Any}'><span class="jlbinding">BiodiversityObservationNetworks.tilt</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tilt(layer, α)
```


Performs logistic-exponential tilting on the on a layer with scaling factor α. This is useful for adjusting inclusion probabilities to scale more toward (α &gt; 1) or away from (α &lt; 1) extreme values. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/utilities/tilting.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.unique_permutations-Union{Tuple{T}, Tuple{T, Any}} where T' href='#BiodiversityObservationNetworks.unique_permutations-Union{Tuple{T}, Tuple{T, Any}} where T'><span class="jlbinding">BiodiversityObservationNetworks.unique_permutations</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
unique_permutations(x::T, prefix = T()) where {T}
```


Generate all unique permutations for a multiset `x` without repetition of duplicates. Based on StackOverflow (`https://stackoverflow.com/questions/65051953/julia-generate-all-non-repeating-permutations-in-set-with-duplicates`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/cubesampling.jl#L290-L296" target="_blank" rel="noreferrer">source</a></Badge>

</details>

