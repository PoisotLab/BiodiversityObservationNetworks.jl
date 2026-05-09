# # Evaluating BON Design

# BiodiversityObservationNetworks.jl contains several utilities for measuring properties of a BON.
# These properties can be split into two categories: (1) properties that measure _spatial balance_ 
# and (2) properties that measure _environmental balance_. 

# Let's load the packages and a couple functions from StatsBase

using BiodiversityObservationNetworks
using BiodiversityObservationNetworks.StatsBase: mean, std
using CairoMakie
using SpeciesDistributionToolkit
using Random #hide
Random.seed!(123) #hide

# For these examples, we'll use a `30 x 30` matrix of zeros as our domain, but any domain will work here

domain = zeros(30, 30)

# Let's compare spatial balance (using `MoransI`) for three different methods: (1) Simple Random Sampling,
# (2) BalancedAcceptance and (3) GRTS.

# We'll use 500 random BON samples per sampling method
reps = 500

# And use in-line loops to record Moran's I across 500 random BONs in a `50 x 50` raster
grts_I = [spatialbalance(MoransI(), domain, sample(GRTS(), domain)) for r in 1:reps]
bas_I = [spatialbalance(MoransI(), domain, sample(BalancedAcceptance(), domain)) for r in 1:reps]
srs_I =  [spatialbalance(MoransI(), domain, sample(SimpleRandom(), domain)) for r in 1:reps]


#
# fig-morans-i-comparison
f = Figure(size=(750,500))
NUM_BINS = 25
ax = Axis(f[1,1], xlabel = "Moran's I", ylabel = "Number of Samplers")
hist!(ax, srs_I, color=(:purple, 0.7), bins=NUM_BINS, label = "Simple\nRandom")
hist!(ax, grts_I, color=(:dodgerblue, 0.7), bins=NUM_BINS, label = "GRTS")
hist!(ax, bas_I, color=(:green, 0.7), bins=NUM_BINS, label = "Balanced\nAcceptance")
annotation!(ax, 200, 0, -0.3, 60,
    text = "More Balanced",
    style = Ann.Styles.LineArrow()
)
axislegend(position = :rt)
current_figure() #hide


# We'll now compute `VoronoiVariance` in the same way:

grts_vv = [spatialbalance(VoronoiVariance(), domain, sample(GRTS(), domain)) for r in 1:reps]
bas_vv = [spatialbalance(VoronoiVariance(), domain, sample(BalancedAcceptance(), domain)) for r in 1:reps]
srs_vv =  [spatialbalance(VoronoiVariance(), domain, sample(SimpleRandom(), domain)) for r in 1:reps]

# and visualize

#
# fig-vv-comparison
f = Figure(size=(750,500))
NUM_BINS = 25
ax = Axis(f[1,1], xlabel = "Voronoi Variance", ylabel = "Number of Samplers")
hist!(ax, srs_vv, color=(:purple, 0.7), bins=NUM_BINS, label = "Simple\nRandom")
hist!(ax, grts_vv, color=(:dodgerblue, 0.7), bins=NUM_BINS, label = "GRTS")
hist!(ax, bas_vv, color=(:green, 0.7), bins=NUM_BINS, label = "Balanced\nAcceptance")
annotation!(ax, 150, 0, 0.6, 40,
    text = "More Balanced",
    style = Ann.Styles.LineArrow()
)
xlims!(ax, 0,1)
axislegend(position = :rt)
current_figure() #hide

