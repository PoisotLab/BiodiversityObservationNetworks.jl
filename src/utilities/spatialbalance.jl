abstract type SpatialBalanceMetric end 

"""
    MoransI
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

function _inclusion_indicator(raster, bon)
    selected_cart_idx = [CartesianIndex(SDT.SimpleSDMLayers.__get_grid_coordinate_by_crs(raster, node.coordinate...)) for node in bon]
    mat = zeros(Bool, size(raster))
    mat[selected_cart_idx] .= 1
    indic = raster.indices .& mat
    return [indic[idx] for idx in findall(raster.indices)]
end

function spatialbalance(::MoransI, raster::SDMLayer, bon::BiodiversityObservationNetwork)
    a = _inclusion_indicator(raster, bon) 

    Es, Ns = SDT.eastings(raster), SDT.northings(raster)
    coords = hcat([[j,i] for i in Ns, j in Es][raster.indices]...)
    
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


spatialbalance(::VoronoiVariance, bon::BiodiversityObservationNetwork, geom) = spatialbalance(VoronoiVariance, bon, geom)

function spatialbalance(
    ::Type{VoronoiVariance}, 
    bon::BiodiversityObservationNetwork,
    geom
)
    vor = voronoi(bon, geom)

    total_inclusion_prob = 50 .* GO.area.(vor) ./ GO.area(geom)
    return sum((total_inclusion_prob .- 1).^2)
end
