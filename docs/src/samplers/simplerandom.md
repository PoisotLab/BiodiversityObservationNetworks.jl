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