module BiodiversityObservationNetworks
    using Clustering
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
    using SpecialFunctions
    using SpeciesDistributionToolkit
    using StatsBase 
    using Random

    import GeoInterface as GI
    import GeometryOps as GO
    import SpeciesDistributionToolkit as SDT
    import SpeciesDistributionToolkit.GeoJSON as GJSON
    import SpeciesDistributionToolkit.SimpleSDMLayers.ArchGDAL as AGDAL
    import DelaunayTriangulation as DT

    export BiodiversityObservationNetwork
    export Node
    export Polygon
    export Raster
    export RasterStack

    export BONSampler
    export MultistageSampler
    export SimpleRandom
    export Grid
    export KMeans
    export CubeSampling
    export SpatiallyStratified
    export BalancedAcceptance
    export GeneralizedRandomTessellatedStratified
    export AdaptiveHotspot
    export UncertaintySampling
    export Pivotal

    export sample
    export datatype
    export nonempty
    export is_polygonizable, is_rasterizable, is_bonifyable

    export voronoi
    
    export cornerplot, bonplot

    include(joinpath("geometry", "bon.jl"))
    include(joinpath("geometry", "polygon.jl"))
    include(joinpath("geometry", "raster.jl"))
    include(joinpath("geometry", "stack.jl"))

    include("sample.jl")
    
    include(joinpath("samplers", "simplerandom.jl"))
    include(joinpath("samplers", "grid.jl"))
    include(joinpath("samplers", "kmeans.jl"))
    include(joinpath("samplers", "cube.jl"))
    include(joinpath("samplers", "spatiallystratified.jl"))
    include(joinpath("samplers", "balancedacceptance.jl"))
    include(joinpath("samplers", "grts.jl"))
    include(joinpath("samplers", "pivotal.jl"))

    include(joinpath("samplers", "adaptivehotspot.jl"))
    include(joinpath("samplers", "uncertainty.jl"))

    include(joinpath("utilities", "voronoi.jl"))

    include("plotting.jl")


end     
