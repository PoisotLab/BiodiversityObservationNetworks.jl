module BiodiversityObservationNetworks
    using SimpleSDMLayers
    using Distributions
    using NeutralLandscapes 
    using Random
    using HaltonSequences
    using StatsBase
    using SpecialFunctions
    using Statistics
    using Base: @kwdef

    include("types.jl")
    export SpatialSampler

    include("sampler.jl")
    export rand, rand!

    include("balancedacceptance.jl")
    export BalancedAcceptance

    include("adaptivespatialsampling.jl")
    export AdaptiveSpatialSampling

    include("_helpers.jl")
    export makesdm, makeoccurrence, makebon, makeenv

end

# makebon(makesdm())

