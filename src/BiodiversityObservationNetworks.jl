module BiodiversityObservationNetworks

using Distributions
using Random
using HaltonSequences
using StatsBase
using SpecialFunctions
using ProgressMeter
using SliceMap

include("types.jl")
export BONSeeder, BONRefiner, BONSampler

include("balancedacceptance.jl")
export BalancedAcceptance

include("adaptivespatial.jl")
export AdaptiveSpatial

include("uniqueness.jl")
export Uniqueness


include("seed.jl")
export seed, seed!

include("refine.jl")
export refine, refine!

include("entropize.jl")
export entropize, entropize!

include("optimize.jl")
export optimize

include("utils.jl")
export stack, squish, _squish, _norm

using Requires
function __init__()
    @require SpeciesDistributionToolkit="72b53823-5c0b-4575-ad0e-8e97227ad13b" include(joinpath("integrations", "simplesdms.jl"))
    @require Zygote="e88e6eb3-aa80-5325-afca-941959d7151f" include(joinpath("integrations", "zygote.jl"))
end




end
