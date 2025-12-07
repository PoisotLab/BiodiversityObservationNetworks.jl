# # Getting Started with BiodiversityObservationNetworks.jl

using BiodiversityObservationNetworks
import Random
using CairoMakie
Random.seed!(1234567890); #hide
CairoMakie.activate!(; px_per_unit = 3) #hide

# In this tutorial, we will cover the basics of how to use BiodiversityObservationNetworks.jl. We'll start by loading the package.

# ## _Hello World_ in BiodiversityObservationNetworks.jl

# The primary use of BiodiversityObservationNetworks is for generating [`BiodiversityObservationNetwork`](@ref) using a variety of point-selection algorithms. Generating networks, regardless of the specific algorithm chosen, is done using the [`sample`](@ref) method. The simplest point-selection algorithm is [`SimpleRandom`](@ref), where each location in space (what is meant be location? More on this in the next section) has an equal probability of inclusion.

bon = sample(SimpleRandom())

# By default, [`SimpleRandom`](@ref) (and every other sampler), chooses 50 points. Without any other arguments, [`sample`](@ref) chooses points from a raster than covers the entire globe, with each pixel representing a 1˚ by 1˚ region.

# When a version of the [Makie](https://docs.makie.org/v0.22/) package for data visualization is loaded, we can use built-in functions to visualize the network. We'll use the `CairoMakie` backend. (Learn mroe about Makie backends [here](https://docs.makie.org/stable/explanations/backends/backends#What-is-a-backend))

# We can adjust the number of points to generate by passing an integer directly to [`SimpleRandom`](@ref), i.e.

srs = SimpleRandom(150)

# Alternatively, we can use the the `num_nodes` keyword argument

srs = SimpleRandom(num_nodes=150)

# Both of these methods for adjusting the number of nodes is supported for all sampling algorithms.
# Let's sample and visualize a [`BiodiversityObservationNetwork`](@ref) with more points
# One thing you may notice about the [`BiodiversityObservationNetwork`](@ref) generated using [`SimpleRandom`](@ref) is that many of the points are clumped together. Many of the sampling algorithms in BiodiversityObservationNetworks aim to produce points that are _spatially balanced_, meaning they are well spread out across space, with little clumping. 
# One such sampler is [`BalancedAcceptance`](@ref). Let's similarly make a [`BiodiversityObservationNetwork`](@ref) with 300 nodes that are spatially balanced.

bon = sample(BalancedAcceptance(150))

# Much better! However, we are still missing some crucial things here. For example, what if we only want to select sites on land? This brings us to applying sampling algorithms to different _geometries_.

# ## Geometries in BiodiversityObservationNetworks.jl

# For most practical use-cases, we aren't interested in developing a [`BiodiversityObservationNetwork`](@ref) for all of Earth, but instead for a small region delineated by a polygon, or represented using raster data, which may contain useful covariate information that we want to incorporate into our BON design.  
# In this case, we want our sampling algorithm `algo` to work on some `geometry` which specifies the spatial domain from which sites should be selected. In this case, we still use `sample`, and pass the spatial domain `geometry` as the second argument, i.e. `sample(algo, geometry)`


# ## A Polygon as a Geometry
# For example, maybe we only want to draw points on land. consider drawing a spatially balanced sample using [`BalancedAcceptance`](@ref) for the nation of Colombia. We can start by downloading a Polygon representing the land

using SimpleSDMPolygons

# Then, we can download a polygon from the NaturalEarth database that represents the land on Earth .

land = getpolygon(PolygonData(NaturalEarth, Land))

# and we can plot it to confirm it's what we expect


# fig-foo
lines(land)
current_figure() #hide

# Now we can generate a [`BiodiversityObservationNetwork`](@ref) using [`BalancedAcceptance`](@ref) in the same way as before, but while passing `col` as a second argument to [`sample`](@ref)

bon = sample(BalancedAcceptance(), land)

# Wahoo 🥳. We've done it. 

# ## A Raster as a Geometry

# Okay, but what if you've got raster data that describes useful environmental covariates? Or a mask of where we can sample? We can use that too.

# Let's start by downloading a raster to use as a source. 


# ## Can we feed a BON to itself?

# Ethically I'm not 100% sure. But it is technically possible. That's both true about sampling BONs from BONs, and the moral of Jurassic Park (1994). 
# Let's download Switzerland.

swi = getpolygon(PolygonData(OpenStreetMap, Places), place="Switzerland")

# and now lets choose a buncha random places in there

candidate_bon = sample(SimpleRandom(500), swi)

# Wow. We're doing groundbreaking work here.
# Next up, let's choose a set of spatially balanced coordinates from this set of candidates. We'll do this using a different sampling algorithm, called the Pivotal method [Grafstrom2012SpaBal](@cite), [`Pivotal`](@ref). Is this because [`BalancedAcceptance`](@ref) doesn't work on point-like geometries? Yes

num_points_to_pick = 30
sample(Pivotal(num_points_to_pick), candidate_bon)

