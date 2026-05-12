module BiodiversityObservationNetworks
    using Random 
    using StatsBase
    using TestItems
    using Distributions
    using LinearAlgebra
    using NearestNeighbors
    using SpecialFunctions
    using Distances
    using HaltonSequences
    using JuMP
    using HiGHS
    using Crayons
    using MultivariateStats
    using DelaunayTriangulation

    if Sys.WORD_SIZE == 64
        const _FLOAT_TYPE = Float64
    else
        const _FLOAT_TYPE = Float32
    end

    include("sampler.jl")
    include("bon.jl")
    include("pool.jl")
    include("sample.jl")

    include(joinpath("utilities", "nearestneighbors.jl"))

    include(joinpath("samplers", "simplerandom.jl"))
    include(joinpath("samplers", "spatiallycorrelatedpoisson.jl"))
    include(joinpath("samplers", "pivotal.jl"))
    include(joinpath("samplers", "balancedacceptance.jl"))
    include(joinpath("samplers", "grts.jl"))
    include(joinpath("samplers", "adaptivehotspot.jl"))
    include(joinpath("samplers", "cubesampling.jl"))
    include(joinpath("samplers", "stratified.jl"))

    include(joinpath("utilities", "rarity.jl"))
    include(joinpath("utilities", "velocity.jl"))
    include(joinpath("utilities", "evaluation.jl"))

    include("show.jl")

    export sample 
    export CandidatePool
    export BiodiversityObservationNetwork

    export BONSampler
    export SimpleRandom, SpatiallyCorrelatedPoisson, Pivotal, BalancedAcceptance, GRTS, AdaptiveHotspot, CubeSampling, Stratified

    export rarity
    export RarityMetric, DistanceToMedian, MultivariateEnvironmentalSimilarity, DistanceToAnalogNode, WithinRange

    export velocity
    export ClosestAnalogue, Loarie2009

    export spatialbalance, evaluate
    export SamplingMetric, MoransI, VoronoiVariance, JensenShannon
end     
