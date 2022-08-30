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
export squish, optimize

using Requires
function __init__()
    @require NeutralLandscapes="71847384-8354-4223-ac08-659a5128069f" include(joinpath("integrations", "neutrallandscapes.jl"))
    @require SimpleSDMLayers="2c645270-77db-11e9-22c3-0f302a89c64c" include(joinpath("integrations", "simplesdms.jl"))
    @require Zygote="e88e6eb3-aa80-5325-afca-941959d7151f" include(joinpath("integrations", "zygote.jl"))
    @require Optim="429524aa-4258-5aef-a3af-852621145aeb" include(joinpath("integrations", "optim.jl"))
end




end
