# # Targeting Unique and High-Velocity Climates for Sampling

# This tutorial will focus on using the utilities in BiodiversityObservationNetworks.jl to use environmental data to compute rarity (how unique the environmental conditions in a location are) and velocity (how quickly environmental conditions are changing), and using those results to target spatially balanced samples toward rare climates, high-velocity climates, or both. 

# ## Acquiring climate data

# We will use SpeciesDistributionToolkit.jl to download climate data. 

using Pkg

@info Pkg.status()

using BiodiversityObservationNetworks
using SpeciesDistributionToolkit
using CairoMakie

# We'll start by downloading a polygon for the state of Washington state, which is the region we will use for this tutorial.

aoi = getpolygon(PolygonData(OpenStreetMap, Places), place = "Washington State")

# We can then load the 19 Bioclimatic variables from CHELSA as follows. 

bioclim = [SDMLayer(RasterData(CHELSA1, BioClim); SpeciesDistributionToolkit.SimpleSDMPolygons.boundingbox(aoi)..., layer = i) for i in 1:19]

# Further documentation on the various available data can be found in [SpeciesDistributionToolkit.jl's documentation](https://poisotlab.github.io/SpeciesDistributionToolkit.jl/v1.9.1/manual/retrieval/list-raster-layers).

# We'll now mask the bioclimatic layers to only include regions inside Washington:

mask!(bioclim, aoi)

# And now we'll visualize the first layer, which is the mean annual temperature:

#
# fig-washington-temp
f = Figure() 
ax = Axis(f[1,1], aspect=DataAspect())
heatmap!(ax, bioclim[1])
lines!(ax, aoi)
current_figure() #hide


# ## Measuring Climate Rarity

# BiodiversityObservationNetworks.jl contains several utilities for quantifying how rare the environmental conditions at a particular location are.

# ### Distance to median climate

# The simplest is [`DistanceToMedian`](@ref), which is each pixel's distance in environmental space to the median environmental condiations.
rar = evaluate(
    DistanceToMedian(), 
    bioclim
)
heatmap(rar)



# ### MultivariateEnvironmentalSimilarity

# Consider 

mess = evaluate(
    MultivariateEnvironmentalSimilarity(),
    bioclim
)
heatmap(quantize(mess))



# ### Distance to Analog Node

# Physical distance to closest node in environment space

bon = sample(SimpleRandom(), bioclim)

rar = evaluate(
    DistanceToAnalogNode(), 
    bon,
    bioclim;
)

heatmap(quantize(rar))
scatter!(bon, color=:red)
current_figure()


rar = evaluate(
    WithinRange(), 
    bon,
    bioclim;
)
heatmap(rar)


# ## Targeting Rare Climates for Sampling 

rar = evaluate(
    DistanceToMedian(), 
    bioclim
)

# # Target rare environments with BAS w/ inclusion 
αs = [5, 3, 1]
bons = [sample(BalancedAcceptance(), rar, inclusion=exp.(αs[i] * rar)) for i in eachindex(αs)]

#
# fig-bas-rarity
f = Figure(size=(1000, 400))
axes = [Axis(f[1,i], title = "α = $(αs[i])", aspect=DataAspect()) for i in 1:3]
for i in eachindex(αs)
    heatmap!(axes[i], rar)
    scatter!(axes[i], bons[i], color=:white, strokewidth=1, strokecolor=:black)
end
current_figure() #hide


# ## Velocity

# Loarie 2009 climate velocity cite.

future_bioclim = [
    SDMLayer(RasterData(CHELSA1, BioClim), Projection(RCP45, ACCESS1_0); 
    SpeciesDistributionToolkit.SimpleSDMPolygons.boundingbox(aoi)..., layer = i) 
    for i in 1:19
]
mask!(future_bioclim, aoi)


vel = evaluate(Loarie2009(), [2000, 2050], [bioclim[1], future_bioclim[1]])
heatmap(vel)

vel = evaluate(Loarie2009(), [2000, 2050], [bioclim, future_bioclim])
heatmap(vel)

