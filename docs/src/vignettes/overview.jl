# # An introduction to BiodiversityObservationNetworks

# In this vignette, we will walk through the basic functionalities of the
# package, by generating a random uncertainty matrix, and then using a *seeder*
# and a *refiner* to decide which locations should be sampled in order to gain
# more insights about the process generating this entropy.

using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots

# In order to simplify the process, we will use the *NeutralLandscapes* package
# to generate a 100×100 pixels landscape, where each cell represents the entropy
# (or information content) in a unit we can sample:

U = rand(MidpointDisplacement(0.5), (100, 100))
heatmap(U'; aspectratio = 1, frame = :none, c = :lapaz)

# In practice, this uncertainty matrix is likely to be derived from an
# application of the hyper-parameters optimization step, which is detailed in
# other vignettes.

# The first step of defining a series of locations to sample is to use a
# `BONSeeder`, which will generate a number of relatively coarse proposals that
# cover the entire landscape, and have a balanced distribution in space. We do
# so using the `BalancedAcceptance` sampler, which can be tweaked to capture
# more (or less) uncertainty. To start with, we will extract 200 candidate
# points, *i.e.* 200 possible locations which will then be refined.

pack = seed(BalancedAcceptance(; numpoints = 200), U);

# The output of a `BONSampler` (whether at the seeding or refinement step) is
# always a tuple, storing in the first position a vector of `CartesianIndex`
# elements, and in the second position the matrix given as input. We can have a
# look at the first five points:

first(pack)[1:5]

# Although returning the input matrix may seem redundant, it actually allows to
# chain samplers together to build pipelines that take a matrix as input, and
# return a set of places to sample as outputs; an example is given below.

# The positions of locations to sample are given as a vector of
# `CartesianIndex`, which are coordinates in the uncertainty matrix. Once we
# have generated a candidate proposal, we can further refine it using a
# `BONRefiner` -- in this case, `AdaptiveSpatial`, which performs adaptive
# spatial sampling (maximizing the distribution of entropy while minimizing
# spatial auto-correlation).

candidates, uncertainty = pack
locations, _ = refine(candidates, AdaptiveSpatial(; numpoints = 50), uncertainty)
locations[1:5]

# The reason we start from a candidate set of points is that some algorithms
# struggle with full landscapes, and work much better with a sub-sample of them.
# There is no hard rule (or no heuristic) to get a sense for how many points
# should be generated at the seeding step, and so experimentation is a must!

# The previous code examples used a version of the `seed` and `refine` functions
# that is very useful if you want to change arguments between steps, or examine
# the content of the candidate pool of points. In addition to this syntax, both
# functions have a curried version that allows chaining them together using
# pipes (`|>`):

locations =
    U |>
    seed(BalancedAcceptance(; numpoints = 200)) |>
    refine(AdaptiveSpatial(; numpoints = 50)) |>
    first

locations[1:5]

# This works because `seed` and `refine` have curried versions that can be used
# directly in a pipeline. Proposed sampling locations can then be overlayed onto the original uncertainty matrix:

plt = heatmap(U'; aspectratio = 1, frame = :none, c = :lapaz)
scatter!(plt, [x[1] for x in locations], [x[2] for x in locations], ms=2.5, mc=:white, label="")