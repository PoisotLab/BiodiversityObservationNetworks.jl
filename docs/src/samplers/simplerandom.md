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
using GeoMakie
import SpeciesDistributionToolkit as SDT
```

now sample a [`BiodiversityObservationNetwork`](@ref)

```@example 1
num_nodes = 50
corsica = gadm("FRA", "Corse")
bon = sample(SimpleRandom(num_nodes), corsica)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
f
```

### Using a vector of Polygons

First, load France and its administrative regions:

```@example 1
france_states = gadm("FRA", 1)
```

Now sample a BON. For a vector of polygons, [`SimpleRandom`](@ref) will sample
from each polygon separately. Here we'll sample 10 nodes in each administrative region.


```@example 1
bon = sample(SimpleRandom(10), france_states)
```

And now plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, france_states, axistype=GeoAxis)
f
```

### Using a raster 

First, load mean annual temperature for France using SDT

```@example 1
france = gadm("FRA")
temp = SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=1, SDT.boundingbox(france)...)
temp = SDT.mask!(temp, france)
```

now sample a BON

```@example 1
bon = sample(SimpleRandom(50), temp)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, temp, axistype=GeoAxis)
f
```

### Using a set of rasters 

Sampling with a set of rasters works almost identically:

```@example 1
france = gadm("FRA")
layers = SDT.mask!([SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=i, SDT.boundingbox(france)...) for i in [1,12]], france)
```

now sample a BON

```@example 1
bon = sample(SimpleRandom(50), layers)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, temp, axistype=GeoAxis)
f
```

### Using a BON

[`SimpleRandom`](@ref) also works with [`BiodiversityObservationNetwork`](@ref)
as a domain.

Lets generate 250 points using the Corsica [`Polygon`](@ref) first, and plot for good measure

```@example 1
bon_candidate = sample(SimpleRandom(250), corsica)

f = Figure(size=(500, 500))
bonplot(f[1,1], bon_candidate, corsica, axistype=GeoAxis)
f
```

and then sample a subset of 100 of them

```@example 1
bon = sample(SimpleRandom(100), bon_candidate)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
f
```