module BiodiversityObservationNetworks

using Distributions
using Random
using HaltonSequences
using StatsBase
using SpecialFunctions
using ProgressMeter
using SliceMap
using Base: @kwdef
using Term
using Term.TermMarkdown
using Markdown

include(joinpath("types", "samplers.jl"))
include(joinpath("types", "groups.jl"))
include(joinpath("types", "targets.jl"))
include(joinpath("types", "layers.jl"))
include(joinpath("types", "weights.jl"))

export BONSeeder, BONRefiner, BONSampler
export Layer, LayerSet, Target, Group, Weights
export name
export getlayers, gettargets, getgroups
export numlayers, numtargets, numgroups
   
    
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
    @require SimpleSDMLayers="2c645270-77db-11e9-22c3-0f302a89c64c" include(joinpath("integrations", "simplesdms.jl"))
    @require Zygote="e88e6eb3-aa80-5325-afca-941959d7151f" include(joinpath("integrations", "zygote.jl"))
end




end


