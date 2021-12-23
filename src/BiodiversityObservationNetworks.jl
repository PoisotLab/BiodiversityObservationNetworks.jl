module BiodiversityObservationNetworks
    using SimpleSDMLayers
    using Distributions
    using NeutralLandscapes 
    using Random
    using HaltonSequences
    using Base: @kwdef

    include("types.jl")
    export SpatialSampler

    include("sampler.jl")
    export rand, rand!

    include("balancedacceptance.jl")
    export BalancedAcceptance

    include("_helpers.jl")
    export makesdm, makeoccurrence, makebon, makeenv

end

# makebon(makesdm())

