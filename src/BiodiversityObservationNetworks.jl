module BiodiversityObservationNetworks
    using DelaunayTriangulation
    using Distances
    using Distributions
    using Extents
    using GeometryOps
    using HaltonSequences
    using HiGHS
    using JuMP
    using LinearAlgebra
    using MultivariateStats
    using NearestNeighbors
    using SparseArrays
    using SpecialFunctions
    using SpeciesDistributionToolkit
    using StatsBase
    using TestItems

    const DT = DelaunayTriangulation
    const SB = StatsBase
    const GO = GeometryOps

    abstract type BONSampler end
    allows_custom_inclusion(::BONSampler) = false

    _DEFAULT_NUM_NODES = 50

    export sample
    export RasterDomain, PolygonDomain, RasterStack
    export extent, contains

    export voronoi

    export spatialbalance
    export VoronoiVariance, MoransI

    export rarity
    export DistanceToAnalogNode, WithinRange, MultivariateEnvironmentalSimilarity, DistanceToMedian

    export jensenshannon

    export BiodiversityObservationNetwork

    export BONSampler
    export AdaptiveHotspot, BalancedAcceptance, CubeSampling, GeneralizedRandomTesselated, Gridded, Pivotal, SimpleRandom, SpatiallyCorrelatedPoisson, UncertaintySampling, SpatiallyStratified

    include(joinpath("domains", "raster.jl"))
    include(joinpath("domains", "stack.jl"))
    include(joinpath("domains", "bon.jl"))
    include(joinpath("domains", "polygon.jl"))
    include(joinpath("domains", "conversion.jl"))

    include(joinpath("mask.jl"))
    include(joinpath("inclusion.jl"))

    include(joinpath("utilities", "nearestneighbors.jl"))
    include(joinpath("utilities", "voronoi.jl"))
    include(joinpath("utilities", "spatialbalance.jl"))
    include(joinpath("utilities", "tilting.jl"))
    include(joinpath("utilities", "distances.jl"))
    include(joinpath("utilities", "clustering.jl"))
    include(joinpath("utilities", "rarity.jl"))


    include(joinpath("sample.jl"))
    include(joinpath("samplers", "simplerandom.jl"))
    include(joinpath("samplers", "spatiallycorrelatedpoisson.jl"))
    include(joinpath("samplers", "cubesampling.jl"))
    include(joinpath("samplers", "balancedacceptance.jl"))
    include(joinpath("samplers", "grts.jl"))
    include(joinpath("samplers", "pivotal.jl"))
    include(joinpath("samplers", "adaptivehotspot.jl"))
    include(joinpath("samplers", "stratified.jl"))


end     
