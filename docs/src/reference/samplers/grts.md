# Generalized Random Tessellated Stratified Sampling

```@docs; canonical=false
GeneralizedRandomTessellatedStratified
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
corsica = openstreetmap("Corse")
bon = sample(GeneralizedRandomTessellatedStratified(num_nodes), corsica)
```

and plot

```@example 1
f = Figure(size=(500, 500))
bonplot(f[1,1], bon, corsica, axistype=GeoAxis)
f
```