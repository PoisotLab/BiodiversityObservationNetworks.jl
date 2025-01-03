module BiodiversityObservationNetworks
    using Clustering
    using DelaunayTriangulation
    using Distributions
    using GeometryOps
    using GeoInterface
    using MultivariateStats
    using SpeciesDistributionToolkit
    using StatsBase 
    using Random

    import GeoInterface as GI
    import GeometryOps as GO
    import SpeciesDistributionToolkit as SDT
    import SpeciesDistributionToolkit.GeoJSON as GJSON

    export BiodiversityObservationNetwork
    export Node
    export Polygon
    export Raster
    export RasterStack

    export BONSampler
    export SimpleRandom, Grid, KMeans, SpatiallyStratified
    export sample
    export datatype
    export nonempty
    export is_polygonizable, is_rasterizable, is_bonifyable

    export cornerplot, bonplot

    include(joinpath("geometry", "bon.jl"))
    include(joinpath("geometry", "polygon.jl"))
    include(joinpath("geometry", "raster.jl"))
    include(joinpath("geometry", "stack.jl"))

    include("sample.jl")
    include(joinpath("samplers", "simplerandom.jl"))
    include(joinpath("samplers", "grid.jl"))
    include(joinpath("samplers", "kmeans.jl"))
    include(joinpath("samplers", "spatiallystratified.jl"))

    include("plotting.jl")


end     