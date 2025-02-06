# CubeSampling

```@docs; canonical=false
CubeSampling
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

Load a [`RasterStack`](@ref) of predictors 

```@example 1
france = SDT.gadm("FRA")
layers = RasterStack(SDT.SDMLayer[SDT.mask!(SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=i, SDT.boundingbox(france)...), france) for i in [1,12]])
```


```@example 1
num_nodes = 50
bon = sample(CubeSampling(num_nodes), layers)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, france, axistype=GeoAxis)
f
```