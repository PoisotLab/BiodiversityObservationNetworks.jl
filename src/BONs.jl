module BONs
    using SimpleSDMLayers
    using Distributions
    using NeutralLandscapes 

    include("types.jl")

    include("_helpers.jl")
    export makesdm, makeoccurrence, makebon, makeenv

end

