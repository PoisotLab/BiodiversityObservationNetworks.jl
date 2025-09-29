abstract type SpatialBalanceMetric end 

"""
    spatialbalance
"""
spatialbalance(::Type{T}, bon, domain) where T = spatialbalance(T(), bon, domain)


spatialbalance(
    metric::SpatialBalanceMetric,
    bon::BiodiversityObservationNetwork,
    domain
) = spatialbalance(metric, bon, to_domain(domain))

"""
    MoransI

`MoransI` is a measure of spatial balance proposed by [Tille2018MeaSpa](@cite).

Conceptually, the idea is to use [Moran's I](https://en.wikipedia.org/wiki/Moran%27s_I), a measure of spatial autocorrelation, on an indicator variable ``\\delta_i``, which is ``1`` if unit ``i`` is included in the sample, and 0 otherwise. 

If ``\\delta`` has a negative value, this means included samples have negative autocorrelation and therefore are more spread out than if generated at random. 

In principle, Moran's I should exist on the interval ``[-1,1]``, but the original definition uses a renormalization that doesn't always guarantee this, so [Tille2018MeaSpa](@cite) introduce a slight variation, called ``I_C`` so this always holds.

Note that the performance of this algorithm scales as ``O(N^3)``, where ``N`` is the number of all possible candidate locations, and therefore for large rasters this method is unlikely to be performant. A faster approximation could be made by taking a smaller random sample of the candidate points. 
"""
struct MoransI <: SpatialBalanceMetric end 

function _morans_weight_matrix(N, n, kdtree, coords)
    W = spzeros(N,N)
    k = Int(floor(N / n))
    resid = (N/n) - k
    for i in 1:N
        nn, _ = NearestNeighbors.knn(kdtree, coords[:,i], k+2, true)
        for j in nn[2:end-1]
            W[i,j] = 1
        end 
        W[i,nn[end]] = resid
    end
    return W
end 

function _get_inclusion_indicator(raster::RasterDomain{<:SDMLayer}, bon)
    selected_cart_idx = [CartesianIndex(SpeciesDistributionToolkit.SimpleSDMLayers.__get_grid_coordinate_by_crs(raster.data, node...)) for node in bon]

    mat = zeros(Bool, size(raster))
    mat[selected_cart_idx] .= 1
    indic = raster.pool .& mat
    return [indic[idx] for idx in getpool(raster)]
end

function _get_inclusion_indicator(raster::RasterDomain{<:Matrix}, bon)
    [p in bon.nodes for p in getpool(raster)]
end


"""
    spatialbalance(::MoransI, raster::SDMLayer, bon::BiodiversityObservationNetwork)
"""
function spatialbalance(
    ::MoransI, 
    bon::BiodiversityObservationNetwork,
    raster::RasterDomain
)

    inclusion_indicator = _get_inclusion_indicator(raster, bon) 
    a = inclusion_indicator

    coords = getcoordinates(raster)
    kdtree = NearestNeighbors.KDTree(coords)

    N, n = size(coords,2), size(bon)

    W = _morans_weight_matrix(N, n, kdtree, coords)
    
    _1 = ones(N)

    b1 = N*W'*W 
    b2 = W'*_1 *_1'*W
    B = b1 - b2
    ψ = n/N

    return (a'*W*a - ψ*(_1'*W*a)) / sqrt( ψ*(1-ψ)*a'*B*a )
end 


spatialbalance(
    mi::MoransI, 
    bon::BiodiversityObservationNetwork,
    rs::RasterStack
) = spatialbalance(mi, bon, first(rs.rasters))

"""
    VoronoiVariance

The `VoronoiVariance` method for characterizing the spatial balance of a sample
is based on the initial method proposed by [Stevens2004SpaBal](@cite), and then
extended by [Grafstrom2012SpaCor](@cite).

For a given [`BiodiversityObservationNetwork`](@ref) `bon`, the [Voronoi
tesselation](https://en.wikipedia.org/wiki/Voronoi_diagram) splits the plane
into a series of polygons, where the `i`-th polygon consists of all points in the
plane whose nearest node in `bon` is the `i`-th node.

These polygons can then be used to assess the spatial balance of a sample. 

In an _ideally_ balanced sample, the sum of the inclusion probabilities across
each polygon `i` would equal 1, because in expectation exactly one unit would
be sampled in that region. 

If we define ``v_i`` as the total inclusion probability across all elements of
the population in Voronoi polygon `i`, i.e.

```math
v_i = \\sum_{j \\in i} \\pi_j 
```

then we can assess the spatial balance of a sample by measuring the distance of
``v_i`` from 1 for each polygon.  [Grafstrom2012SpaCor](@cite) proposes the
metric `B`, defined as 

```math
B = \\frac{1}{n} \\sum_{i=1}^n (v_i - 1)^2
```

to measure spatial balance, where *smaller values indicate better spatial balance*.
"""
struct VoronoiVariance <: SpatialBalanceMetric end 

spatialbalance(::VoronoiVariance, bon::BiodiversityObservationNetwork, geom) = spatialbalance(VoronoiVariance, bon, to_polygon(geom))

"""
    spatialbalance(::Type{VoronoiVariance}, bon::BiodiversityObservationNetwork, geom)
"""
function spatialbalance(
    ::VoronoiVariance, 
    bon::BiodiversityObservationNetwork,
    geom::PolygonDomain
)

    vor = voronoi(bon, geom)

    total_inclusion_prob = length(bon) * GO.area.(vor) ./ GO.area(geom.data)
    return sum((total_inclusion_prob .- 1).^2)
end
