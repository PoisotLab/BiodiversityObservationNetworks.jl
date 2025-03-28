# # Package overview

using BiodiversityObservationNetworks
using NeutralLandscapes

# setup a random matrix of uncertainty

U = rand(MidpointDisplacement(0.5), (100, 100))

# seed the initial set of points

candidates, _ = seed(BalancedAcceptance(numpoints=200), U)

# refine the points

locations, _ = refine(candidates, AdaptiveSpatial(numpoints=50), U)

# alternative way using pipes

U |> seed(BalancedAcceptance(numpoints=200)) |> refine(AdaptiveSpatial(numpoints=50)) |> first