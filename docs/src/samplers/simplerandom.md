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

now sample a BON

```@example 1
num_nodes = 50
corsica = SDT.gadm("FRA", "Corse")
bon = sample(SimpleRandom(num_nodes), corsica)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
f
```

### Using a raster 

First, load mean annual temperature for Corsica using SDT

```@example 1
temp = SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=1, SDT.boundingbox(corsica)...)
```

now sample a BON

```@example 1
bon = sample(SimpleRandom(num_nodes), temp)
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
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
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