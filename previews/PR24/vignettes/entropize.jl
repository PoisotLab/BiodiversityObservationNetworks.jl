# # Getting the entropy matrix

# For some applications, we want to place points to capture the maximum amount
# of information, which is to say that we want to sample a balance of *entropy*
# values, as opposed to *absolute* values. In this vignette, we will walk
# through an example using the `entropize` function to convert raw data to
# entropy values.

using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots
using SliceMap

# !!! warning "Entropy is problem-specific"
#     The solution presented in this vignette is a least-assumption solution
#     based on the empirical values given in a matrix of measurements. In a lot
#     of situations, this is not the entropy that you want. For example, if your
#     pixels are storing probabilities of Bernoulli events, you can directly use
#     the entropy of the events in the entropy matrix.

# We start by generating a random matrix of measurements:

measurements = rand(MidpointDisplacement(), (200, 200)) .* 100
heatmap(measurements)

# Using the `entropize` function will convert these values into entropy at the
# pixel scale:

U = entropize(measurements)
heatmap(U)

# The values closest to the median of the distribution have the highest entropy,
# and the values closest to its extrema have an entropy of 0. The entropy matrix
# is guaranteed to have values on the unit interval.

# We can use `entropize` as part of a pipeline, and overlay the points optimized
# based on entropy on the measurement map:

locations =
    measurements |> entropize |> seed(BalancedAcceptance(; numpoints = 500)) |>
    refine(AdaptiveSpatial(; numpoints = 50)) |> first
heatmap(measurements)
scatter!(
    [x[1] for x in locations],
    [x[2] for x in locations];
    ms = 2.5,
    mc = :white,
    label = "",
)