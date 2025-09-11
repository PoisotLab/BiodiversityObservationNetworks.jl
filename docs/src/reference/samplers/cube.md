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

Load predictors

```@example 1
corsica = SDT.getpolygon(PolygonData(OpenStreetMap, Places), place="Corsica")
layers = SDT.mask!([SDT.SDMLayer(SDT.RasterData(SDT.CHELSA2, SDT.BioClim); layer=i, SDT.boundingbox(corsica)...) for i in [1,12]], corsica)
```


```@example 1
num_nodes = 50
bon = sample(CubeSampling(num_nodes), layers)
```
