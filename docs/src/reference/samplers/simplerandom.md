# Simple Random Sampling

```@docs; canonical=false
SimpleRandom
```

## Example 

### Using a polygon 

First, load the packages we will use for this example

```@example 1
using BiodiversityObservationNetworks 
using CairoMakie
using SpeciesDistributionToolkit
import SpeciesDistributionToolkit as SDT
```

now sample a [`BiodiversityObservationNetwork`](@ref)

```@example 1
num_nodes = 50
corsica = SDT.getpolygon(SDT.PolygonData(OpenStreetMap, Places), place="Corsica")
bon = sample(SimpleRandom(num_nodes), corsica)
```

and plot

### Using a raster 

First, load mean annual temperature for Switzerland using SDT

```@example 1
switzerland = SDT.getpolygon(SDT.PolygonData(OpenStreetMap, Places), place="Switzerland")

temp = SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=1, SDT.boundingbox(switzerland)...)
temp = SDT.mask!(temp, switzerland)
```

now sample a BON

```@example 1
bon = sample(SimpleRandom(50), temp)
```

### Using a BON

[`SimpleRandom`](@ref) also works with [`BiodiversityObservationNetwork`](@ref)
as a domain.

Lets generate 250 points using the Corsica polygon first, and plot for good measure

```@example 1
bon_candidate = sample(SimpleRandom(250), corsica)
```

and then sample a subset of 100 of them

```@example 1
bon = sample(SimpleRandom(100), bon_candidate)
```
