module BiodiversityObservationNetworks

using SimpleSDMLayers
using Distributions
using Random
using HaltonSequences
using StatsBase


include("types.jl")
export BONSeeder, BONRefiner, BONSampler

include("balancedacceptance.jl")
export BalancedAcceptance

#include("adaptivespatialsampling.jl")
#export AdaptiveSpatialSampling

include("adaptivespatialsampling.jl")
export AdaptiveSpatialSampling

include("seed.jl")
export seed, seed!

# TODO define these
# export refine, refine!

end
