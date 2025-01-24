# Simple Random Sampling

```@docs; canonical=false
SimpleRandom
```

## Example 

First, load the packages we will use for this example

```@example 
using BiodiversityObservationNetworks 
using CairoMakie
using GeoMakie
import SpeciesDistributionToolkit as SDT
num_nodes = 50
corsica = SDT.gadm("FRA", "Corse")
bon = sample(SimpleRandom(num_nodes), corsica)
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
f
```