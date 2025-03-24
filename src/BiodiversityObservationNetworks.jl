module BiodiversityObservationNetworks
    using Clustering
    using Dates
    using DelaunayTriangulation
    using Distances
    using Distributions
    using GeometryOps
    using GeoInterface
    using HaltonSequences
    using HiGHS
    using JuMP
    using LinearAlgebra
    using MultivariateStats
    using NearestNeighbors
    using SpecialFunctions
    using SpeciesDistributionToolkit
    using StatsBase 
    using TestItems
    using Random

    import GeoInterface as GI
    import GeometryOps as GO
    import SpeciesDistributionToolkit as SDT
    import SpeciesDistributionToolkit.GeoJSON as GJ
    import SpeciesDistributionToolkit.SimpleSDMLayers.ArchGDAL as AG
    import DelaunayTriangulation as DT
    import MultivariateStats as MVStats

    export BiodiversityObservationNetwork
    export Node
    export Polygon
    export Raster
    export RasterStack

    export BONSampler
    export MultistageSampler
    export SimpleRandom
    export Grid
    export CubeSampling
    export SpatiallyStratified
    export BalancedAcceptance
    export WeightedBalancedAcceptance
    export GeneralizedRandomTessellatedStratified
    export AdaptiveHotspot
    export UncertaintySampling
    export SpatiallyCorrelatedPoisson

    export PivotalMethod, KPivotal, KDTreePivotal, ClassicPivotal
    export Pivotal

    export sample
    export datatype
    export nonempty
    export is_polygonizable, is_rasterizable, is_bonifyable

    export cluster
    export KMeans, FuzzyCMeans

    export features
    export jensenshannon
    export voronoi
    
    export velocity
    export VelocityMetric, Loarie2009, ClosestAnalogue

    export rarity
    export RarityMetric, DistanceToMedian, DistanceToAnalogNode, MultivariateEnvironmentalSimilarity

    export spatialbalance
    export MoransI, VoronoiVariance

    export gadm
    export cornerplot, bonplot



    include(joinpath("geometry", "bon.jl"))
    include(joinpath("geometry", "polygon.jl"))
    include(joinpath("geometry", "raster.jl"))
    include(joinpath("geometry", "stack.jl"))
    include(joinpath("geometry", "timeseries.jl"))

    include("sample.jl")
    
    include(joinpath("samplers", "simplerandom.jl"))
    include(joinpath("samplers", "grid.jl"))
    include(joinpath("samplers", "cube.jl"))
    include(joinpath("samplers", "spatiallystratified.jl"))
    include(joinpath("samplers", "balancedacceptance.jl"))
    include(joinpath("samplers", "weightedbas.jl"))
    include(joinpath("samplers", "grts.jl"))
    include(joinpath("samplers", "pivotal.jl"))
    include(joinpath("samplers", "scps.jl"))

    include(joinpath("samplers", "adaptivehotspot.jl"))
    include(joinpath("samplers", "uncertainty.jl"))

    include(joinpath("utilities", "voronoi.jl"))
    include(joinpath("utilities", "clustering.jl"))
    include(joinpath("utilities", "spatialbalance.jl"))
    include(joinpath("utilities", "velocity.jl"))
    include(joinpath("utilities", "rarity.jl"))



    include("overloads.jl")
    include("plotting.jl")
end     
