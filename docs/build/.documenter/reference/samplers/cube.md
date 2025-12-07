
# CubeSampling {#CubeSampling}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.CubeSampling-reference-samplers-cube' href='#BiodiversityObservationNetworks.CubeSampling-reference-samplers-cube'><span class="jlbinding">BiodiversityObservationNetworks.CubeSampling</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CubeSampling
```


`CubeSampling` implements the cube method of (???) for balanced sampling with respect to a set of auxiliary variables (features).

The algorithm proceeds in two phases:
- Flight phase: probabilities are iteratively moved towards 0/1 while preserving linear constraints on the mean of auxiliary variables in expectation (including fixed sample size).
  
- Landing phase: if fractional probabilities remain, an optimization step chooses a 0/1 sample that best matches the target constraints.
  

If `inclusion` is not provided, uniform inclusion probabilities are derived from `sampler.num_nodes` and the domain pool size. Returned nodes correspond to units with final probability 1.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/cf81ee6a79734664b33fd1821e1d682f012a7cc2/src/samplers/cubesampling.jl#L1-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Example {#Example}

### Using a polygon {#Using-a-polygon}

First, load the packages we will use for this example

```julia
using BiodiversityObservationNetworks
using SpeciesDistributionToolkit
using CairoMakie
import SpeciesDistributionToolkit as SDT
```


Load predictors

```julia
corsica = SDT.getpolygon(SDT.PolygonData(OpenStreetMap, Places), place="Corsica")
layers = SDT.mask!([SDT.SDMLayer(SDT.RasterData(SDT.CHELSA2, SDT.BioClim); layer=i, SDT.boundingbox(corsica)...) for i in [1,12]], corsica)
```


```ansi
2-element Vector{SimpleSDMLayers.SDMLayer{UInt16}}:
 🗺️  A 205 × 124 layer (13685 UInt16 cells)
 🗺️  A 205 × 124 layer (13685 UInt16 cells)
```


```julia
num_nodes = 50
bon = sample(CubeSampling(num_nodes), layers)
```


```ansi
BiodiversityObservationNetwork with 50 nodes
```

