
# Balanced Acceptance Sampling {#Balanced-Acceptance-Sampling}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.BalancedAcceptance-reference-samplers-balancedacceptance' href='#BiodiversityObservationNetworks.BalancedAcceptance-reference-samplers-balancedacceptance'><span class="jlbinding">BiodiversityObservationNetworks.BalancedAcceptance</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BalancedAcceptance
```


`BalancedAcceptance` implements Balanced Acceptance Sampling (BAS), which uses Halton sequences to provide spatially well-spread samples. When inclusion probabilities are supplied, a 3D Halton sequence is used: the first two dimensions spread locations spatially and the third acts as a random threshold for acceptance against the inclusion surface.

If inclusion probabilities are not provided, a 2D Halton sequence is used to generate a spatially balanced sample of size `num_nodes` while respecting raster masks.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/cf81ee6a79734664b33fd1821e1d682f012a7cc2/src/samplers/balancedacceptance.jl#L1-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Example {#Example}

### Using a polygon {#Using-a-polygon}

First, load the packages we will use for this example

```julia
using BiodiversityObservationNetworks
using CairoMakie
import SpeciesDistributionToolkit as SDT
```


now sample a [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork)

```@example 1
num_nodes = 50
corsica = SDT.getpolygon(SDT.PolygonData(OpenStreetMap, Places), place="Corsica")
bon = sample(BalancedAcceptance(num_nodes), corsica)
```

