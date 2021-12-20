module BONs
    using SimpleSDMLayers
    using Distributions
    using NeutralLandscapes 
    using Random

    include("types.jl")
    export SpatialSampler, FractalTriad

    include("sampler.jl")
    export rand, rand!

    include("fractaltriad.jl")

    include("_helpers.jl")
    export makesdm, makeoccurrence, makebon, makeenv

end

# makebon(makesdm())

