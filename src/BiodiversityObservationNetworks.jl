module BiodiversityObservationNetworks

using SimpleSDMLayers
using Distributions
using Random
using HaltonSequences
using StatsBase
using SpecialFunctions
using ProgressMeter
using SliceMap
using JuMP
using HiGHS
using LinearAlgebra
using Term
using TestItems

include("types.jl")
export BONSampler
export Sites, numsites, pool
export LayerType, DataLayer, InclusionProbability
export Layer, Stack

include("sample.jl")
export sample

include("exceptions.jl")
export BONException, TooFewSites, TooManySites

include("simplerandom.jl")
export SimpleRandom

include("spatialstratified.jl")
export SpatiallyStratified

include("balancedacceptance.jl")
export BalancedAcceptance

include("weightedbas.jl")
export WeightedBalancedAcceptance

include("adaptivehotspot.jl")
export AdaptiveHotspot

include("cubesampling.jl")
export CubeSampling

include("fractaltriad.jl")
export FractalTriad

include("grts.jl")
export GeneralizedRandomTessellatedStratified

include("uniqueness.jl")
export Uniqueness

#include("seed.jl")
#export seed, seed!

#include("refine.jl")
#export refine, refine!

include("entropize.jl")
export entropize, entropize!

include("utils.jl")
export stack

end
