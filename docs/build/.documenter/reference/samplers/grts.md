
# Generalized Random Tessellated Stratified Sampling {#Generalized-Random-Tessellated-Stratified-Sampling}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.GeneralizedRandomTesselated-reference-samplers-grts' href='#BiodiversityObservationNetworks.GeneralizedRandomTesselated-reference-samplers-grts'><span class="jlbinding">BiodiversityObservationNetworks.GeneralizedRandomTesselated</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/cf81ee6a79734664b33fd1821e1d682f012a7cc2/src/samplers/grts.jl#L1-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Example {#Example}

### Using a polygon {#Using-a-polygon}

First, load the packages we will use for this example

```julia
using BiodiversityObservationNetworks
using CairoMakie
using SpeciesDistributionToolkit
import SpeciesDistributionToolkit as SDT
```


now sample a [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork)

```julia
num_nodes = 50
corsica = SDT.getpolygon(SDT.PolygonData(OpenStreetMap, Places), place="Corsica")
bon = sample(GeneralizedRandomTesselated(num_nodes), corsica)
```


```ansi
BiodiversityObservationNetwork with 50 nodes
```

