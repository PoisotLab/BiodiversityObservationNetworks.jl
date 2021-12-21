module BiodiversityObservationNetworks
    using SimpleSDMLayers
    using Distributions
    using NeutralLandscapes 
    using Random
    using Base: @kwdef

    include("types.jl")
    export SpatialSampler

    include("sampler.jl")
    export rand, rand!

    include("fractaltriad.jl")
    export FractalTriad

    include("_helpers.jl")
    export makesdm, makeoccurrence, makebon, makeenv

end

# makebon(makesdm())

