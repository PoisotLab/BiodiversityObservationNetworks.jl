# Generalized Random Tessellated Stratified Sampling

```@docs; canonical=false
GeneralizedRandomTesselated
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
bon = sample(GeneralizedRandomTesselated(num_nodes), corsica)
```

