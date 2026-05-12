# # Targeting Unique Climates for Sampling

# This tutorial will focus on using the utilities in BiodiversityObservationNetworks.jl to use environmental data to compute rarity (how unique the environmental conditions in a location are).

# It will also demonstrate how inclusion probabilities can be used to target spatially balanced samples toward rare climates, which is applicable more broadly toward targeting any specific variable of interest. 

# ## Acquiring climate data

# We will use SpeciesDistributionToolkit.jl to download climate data. 

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
hidedecorations!(ax)
hidespines!(ax)
heatmap!(ax, bioclim[1])
lines!(ax, aoi)
current_figure() #hide


# ## Measuring Climate Rarity

# BiodiversityObservationNetworks.jl contains several utilities for quantifying how rare the environmental conditions at a particular location are.

# ### Distance to median climate

# The simplest is [`DistanceToMedian`](@ref), which is each pixel's distance in environmental space to the median environmental conditions across the domain.

# All rarity, velocity, and evaluation metrics are built around the [`evaluate`](@ref) method.

# The API is structured similarly to [`sample`](@ref).

# ```julia
# evaluate(::SamplingMetric, domain)
# ```

# or 

# ```julia
# evaluate(::SamplingMetric, domain, bon)
# ```

# for metrics that involve a BON in addition to the domain.

# For [`DistanceToMedian`](@ref), this looks like

rar = evaluate(
    DistanceToMedian(), 
    bioclim
)


# Which we can then visualize

#
# fig-dist-to-med
f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
hm = heatmap!(ax, rar)
Colorbar(f[2,1], hm, vertical = false, label="Distance to Median Environment")
current_figure() #hide


# ### MultivariateEnvironmentalSimilarity

# A second metric for environmental uniqueness is [`MultivariateEnvironmentalSimilarity`], which is derived from [Mesgaran2014HereBe](@cite).

# This is a _similarity_ score, so **higher values indicate similarity, lower values indicate rarity**.

# This can be run via 

mess = evaluate(
    MultivariateEnvironmentalSimilarity(),
    bioclim
)

# and visualized

#
# fig-mess
f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hm = heatmap!(ax, mess)
hidedecorations!(ax)
hidespines!(ax)
Colorbar(f[2,1], hm, vertical = false, label="Similarity to typical environment")
current_figure() #hide


# ### Distance to Analog Node

# Another metric for climate rarity is comparing the environmental conditions across a domain to the conditions at an existing [`BiodiversityObservationNetwork`](@ref).

# [`DistanceToAnalogNode`](@ref) measures the distance of every pixel in the domain to the closest BON site in _environmental space_.

# We can demonstrate this by first sampling a BON

bon = sample(SimpleRandom(), bioclim)

# and running

dist_to_analog = evaluate(
    DistanceToAnalogNode(), 
    bioclim,
    bon
)

# which we can then visualize

# fig-dist-to-analog
f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
hm = heatmap!(ax, dist_to_analog)
Colorbar(f[2,1], hm, vertical = false, label="Environmental distance to closest proxy node")
current_figure() #hide


# ### Within Range

# A simpler method for assessing a BON's coverage in environmental space is [`WithinRange`](@ref), which determines if each pixel in the domain is within the extrema of the environmental variables at each BON -- i.e. if the environmental conditions at a pixel are at least "bounded" by those at the BON sites.

withinrange = evaluate(
    WithinRange(), 
    bioclim,
    bon
)

#
# fig-dist-to-analog
f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hidespines!(ax)
hidedecorations!(ax)
hm = heatmap!(ax, withinrange, colormap=[:grey80, :seagreen4])
scatter!(ax, bon, color=:white, strokewidth=1, strokecolor=:black)
Legend(f[2,1], [PolyElement(color=:grey80), PolyElement(color=:seagreen4)], ["Outside Range", "Within Range"], orientation=:horizontal)
current_figure() #hide



# ## Targeting Rare Climates for Sampling 

# We can use these climate rarity metrics to target regions of high rarity for sampling.

# We will do this using [`BalancedAcceptance`](@ref) with custom inclusion probabilities, which can be used more generally to target any variable of interest.

# We'll start by using [`DistanceToMedian`](@ref) to compute rarity.

rar = evaluate(
    DistanceToMedian(), 
    bioclim
)

# We can sample a BON with inclusion proportional to rarity as follows

bon = sample(BalancedAcceptance(), rar, inclusion = rar)

# Let's visualize this

#
# fig-bas-rarity
f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hidespines!(ax)
hidedecorations!(ax)
heatmap!(ax, rar)
scatter!(ax, bon, color=:white, strokewidth=1, strokecolor=:black)
Colorbar(f[2,1], hm, vertical = false, label="Distance to Median Environment")
current_figure() #hide


# We see the points are generally more concentrated toward the left side of the map, where the environmental conditions are more rare. 

# Pragmatically, this may not be skewed "enough" toward rare regions, so we can adjust the relative inclusion probabilities by applying an exponential transformation with a parameter $\alpha$ adjustsing the "strength" of this transform.

αs = [1,3,5,10]
bons = [sample(BalancedAcceptance(), rar, inclusion=exp.(αs[i] * rar)) for i in eachindex(αs)]

# We can then visualize to see the impact of this transform

#
# fig-bas-alpha-comparison
f = Figure(size=(1000, 600))
axes = [Axis(f[ci[2],ci[1]], title = "α = $(αs[i])", aspect=DataAspect()) for (i,ci) in enumerate(CartesianIndices((1:2,1:2)))]
for i in eachindex(axes)
    heatmap!(axes[i], rar)
    scatter!(axes[i], bons[i], color=:white, strokewidth=1, strokecolor=:black)
end
hidespines!.(axes)
hidedecorations!.(axes)
current_figure() #hide


# ## Next steps

# In the [next tutorial](/tutorials/evaluation), we'll cover how to evaluate the spatial balance and environmental representativeness of a particular BON.