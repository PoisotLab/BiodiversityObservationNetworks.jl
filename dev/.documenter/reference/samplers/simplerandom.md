
# Simple Random Sampling {#Simple-Random-Sampling}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.SimpleRandom-reference-samplers-simplerandom' href='#BiodiversityObservationNetworks.SimpleRandom-reference-samplers-simplerandom'><span class="jlbinding">BiodiversityObservationNetworks.SimpleRandom</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SimpleRandom
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/a7d32bc5f6558ea22e251490da729de69bdb2c0a/src/samplers/simplerandom.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

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
bon = sample(SimpleRandom(num_nodes), corsica)
```


```ansi
BiodiversityObservationNetwork with 50 nodes
```


and plot

### Using a raster {#Using-a-raster}

First, load mean annual temperature for Switzerland using SDT

```julia
switzerland = SDT.getpolygon(SDT.PolygonData(OpenStreetMap, Places), place="Switzerland")

temp = SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=1, SDT.boundingbox(switzerland)...)
temp = SDT.mask!(temp, switzerland)
```


```ansi
🗺️  A 13 × 28 layer with 179 Float32 cells
   Projection: +proj=longlat +datum=WGS84 +no_defs
```


now sample a BON

```julia
bon = sample(SimpleRandom(50), temp)
```


```ansi
BiodiversityObservationNetwork with 50 nodes
```


### Using a BON {#Using-a-BON}

[`SimpleRandom`](/reference/api#BiodiversityObservationNetworks.SimpleRandom) also works with [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork) as a domain.

Lets generate 250 points using the Corsica polygon first, and plot for good measure

```julia
bon_candidate = sample(SimpleRandom(250), corsica)
```


```ansi
BiodiversityObservationNetwork with 250 nodes
```


and then sample a subset of 100 of them

```julia
bon = sample(SimpleRandom(100), bon_candidate)
```


```ansi
BiodiversityObservationNetwork with 100 nodes
```

