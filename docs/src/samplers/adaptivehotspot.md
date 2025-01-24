# Adaptive Hotspot


```@docs; canonical=false
AdaptiveHotspot
```

## Example 

### Using a raster 

First, load the packages we will use for this example

```@example 1
using BiodiversityObservationNetworks 
using CairoMakie
using GeoMakie
using NeutralLandscapes
import SpeciesDistributionToolkit as SDT
```

Now we'll generate a synthetic uncertainty layer using
[`NeutralLandscapes.jl`](http://docs.ecojulia.org/NeutralLandscapes.jl/dev/).


```@example 1
uncertainty = rand(MidpointDisplacement(), (100, 100))
```

Let's plot it

```@example 1
heatmap(uncertainty)
```

and now we sample a BON

```@example 1
bon = sample(AdaptiveHotspot(), uncertainty)
```

and plot it

```@example 1
f = Figure()
bonplot(f[1,1], bon, uncertainty)
f
```

