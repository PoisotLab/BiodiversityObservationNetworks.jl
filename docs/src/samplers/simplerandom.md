# Simple Random Sampling

```@docs; canonical=false
SimpleRandom
```

## Example 

First, load the packages we will use for this example

```@example 1
using BiodiversityObservationNetworks 
using CairoMakie
using GeoMakie
import SpeciesDistributionToolkit as SDT
```

Using a Polygon as the domain

```@example 1
num_nodes = 50
corsica = SDT.gadm("FRA", "Corse")
bon = sample(SimpleRandom(num_nodes), corsica)
```

Now plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
f
```