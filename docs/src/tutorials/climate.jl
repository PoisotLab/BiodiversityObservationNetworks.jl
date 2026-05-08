# # Targeting Unique and High-Velocity Climates for Sampling

# This tutorial will focus on using the utilities in BiodiversityObservationNetworks.jl to use environmental data to compute rarity (how unique the environmental conditions in a location are) and velocity (how quickly environmental conditions are changing), and using those results to target spatially balanced samples toward rare climates, high-velocity climates, or both. 

# ## Acquiring climate data

# We will use SpeciesDistributionToolkit.jl to download climate data. 

using Pkg

using BiodiversityObservationNetworks
using SpeciesDistributionToolkit
using CairoMakie

# We'll start by downloading a polygon for the state of Oregon, which is the region we will use for this tutorial.

aoi = getpolygon(PolygonData(OpenStreetMap, Places), place = "Oregon")

# We can then load the 19 Bioclimatic variables from CHELSA as follows. 

bioclim = [SDMLayer(RasterData(CHELSA1, BioClim); SpeciesDistributionToolkit.SimpleSDMPolygons.boundingbox(aoi)..., layer = i) for i in 1:19]

# Further documentation on the various available data can be found in [SpeciesDistributionToolkit.jl's documentation](https://poisotlab.github.io/SpeciesDistributionToolkit.jl/v1.9.1/manual/retrieval/list-raster-layers).

# We'll now mask the bioclimatic layers to only include regions inside Oregon:

mask!(bioclim, aoi)

# And now we'll visualize the first layer, which is the mean annual temperature:

# fig-oregon-temp
f = Figure() 
ax = Axis(f[1,1], aspect=DataAspect())
heatmap!(ax, bioclim[1])
lines!(ax, aoi)
current_figure() #hide



mess = rarity(
    MultivariateEnvironmentalSimilarity(),
    bioclim
)
heatmap(mess)


rar = rarity(
    DistanceToMedian(), 
    bioclim
)

heatmap(mess)
heatmap(rar)

bon = sample(SimpleRandom(), bioclim)

rar = rarity(
    DistanceToAnalogNode(), 
    bon,
    bioclim;
)

heatmap(rar)
scatter!(bon, color=:red)
current_figure()


rar = rarity(
    WithinRange(), 
    bon,
    bioclim;
)


# ## Velocity

future_bioclim = [
    SDMLayer(RasterData(CHELSA1, BioClim), Projection(RCP45, ACCESS1_0); 
    SpeciesDistributionToolkit.SimpleSDMPolygons.boundingbox(aoi)..., layer = i) 
    for i in 1:19
]
mask!(future_bioclim, aoi)

#vel = velocity(ClosestAnalogue(), bioclim, future_bioclim)


vel = velocity(Loarie2009(), [2000, 2050], [bioclim[12], future_bioclim[12]])

heatmap(vel)
heatmap(quantize(vel))
