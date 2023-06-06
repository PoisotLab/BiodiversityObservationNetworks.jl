# # An introduction to BiodiversityObservationNetworks

# In this vignette, we will walk through the basic functionalities of the
# package, by generating a random uncertainty matrix, and then using a *seeder*
# and a *refiner* to decide which locations should be sampled in order to gain
# more insights about the process generating this entropy.

using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots

# In order to simplify the process, we will use the *NeutralLandscapes* package
# to generate a 100×100 pixels landscape, where each cell represents a unit we
# can sample:

U = rand(MidpointDisplacement(0.5), (100, 100))
heatmap(U; aspectratio = 1, frame = :none, c = :lapaz)

# In practice, this uncertainty matrix is likely to be derived from an
# application of the hyper-parameters optimization step, which is detailed in
# other vignettes.

# The first step of defining a series of locations to sample is to use a
# `BONSeeder`, which will generate a number of relatively coarse proposals that
# cover the entire landscape, and have a balanced distribution in space. We do
# so using the `BalancedAcceptance` sampler, which can be tweaked to capture
# more (or less) uncertainty. To start with, we will extract 200 candidate
# points, *i.e.* 200 possible locations which will then be refined.

pack = seed(BalancedAcceptance(; numpoints = 200), U)
first(pack)[1:5]

# We store the output in `pack`, which is a tuple containing the recommended
# locations and a copy of the uncertainty matrix. This enables to build
# workflows where samplers are chained together.

# The positions of locations to sample are given as a vector of
# `CartesianCoordinates`, which are coordinates of the uncertainty matrix. Once
# we have generated a candidate proposal, we can further refine it using a
# `BONRefiner` -- in this case, `AdaptiveSpatial`, which performs adaptive
# spatial sampling (maximizing the distribution of entropy while minimizing
# spatial auto-correlation).

candidates, uncertainty = pack
locations, _ = refine(candidates, AdaptiveSpatial(; numpoints = 50), uncertainty)

# The reason we start from a candidate set of points is that some algorithms
# struggle with full landscapes, and work much better with a sub-sample of
# them.

# Note that the syntax we would actually use is a lot simpler, and involves
# using pipes (`|>`):

locations =
    U |>
    seed(BalancedAcceptance(; numpoints = 200)) |>
    refine(AdaptiveSpatial(; numpoints = 50)) |>
    first

locations[1:5]

# This works because `seed` and `refine` have curried versions that can be used
# directly in a pipeline.