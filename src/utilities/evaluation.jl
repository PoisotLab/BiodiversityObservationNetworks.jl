"""
    SamplingMetric

Abstract supertype for metrics that evaluate a [`SamplingResult`](@ref).

Implement `evaluate(metric::MyMetric, result::SamplingResult)` to add a
new metric.
"""
abstract type SamplingMetric end


"""
    MoransI <: SamplingMetric

Computes Moran's I on the inclusion indicator variable.

# Description
Calculates spatial autocorrelation of the sample indicator ``\\delta`` (1 if sampled, 0 otherwise).
Negative values indicate spatial inhibition (spread), which is desired for balanced sampling.
"""
struct MoransI <: SamplingMetric end 

function _morans_weight_matrix(N, n, kdtree, coords)
    W = zeros(N,N)
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


function _get_inclusion_indicator(raster::Matrix, bon)
    mat = zeros(Bool, size(raster))
    mat[bon.sites] .= 1
    return [mat[idx] for idx in findall(x-> x isa Number, raster)]
end

function _morans_i(coords, inclusion_indicator, bon)
    a = inclusion_indicator

    kdtree = NearestNeighbors.KDTree(coords)

    N, n = size(coords,2), length(bon)

    W = _morans_weight_matrix(N, n, kdtree, coords)
    
    _1 = ones(N)

    b1 = N*W'*W 
    b2 = W'*_1 *_1'*W
    B = b1 - b2
    ψ = n/N

    return (a'*W*a - ψ*(_1'*W*a)) / sqrt( ψ*(1-ψ)*a'*B*a )
end 



"""
    spatialbalance(::MoransI, domain::Matrix, bon::BiodiversityObservationNetwork)
"""
function spatialbalance(
    ::MoransI, 
    domain::Matrix,
    bon::BiodiversityObservationNetwork,
)
    coords = _FLOAT_TYPE.(hcat([[ci[1], ci[2]] for ci in CartesianIndices(domain)]...))
    inclusion_indicator = _get_inclusion_indicator(domain, bon) 
    _morans_i(coords, inclusion_indicator, bon)
end


"""
    VoronoiVariance <: SamplingMetric

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
struct VoronoiVariance <: SamplingMetric end

"""
    _voronoi_tesselation(bon, domain)

Construct Voronoi polygons for the nodes in a
`BiodiversityObservationNetwork` within the given domain. 
Output polygons are clipped to the domain extent.
"""
function _voronoi_tesselation(bon_coords)
    tri = DelaunayTriangulation.triangulate(bon_coords)
    vor = DelaunayTriangulation.voronoi(tri)

    #polys = [intersect(domain.data, LeanBONsAPI.Polygon(DT.get_polygon_coordinates(vor, i, bbox))) for i in DT.each_polygon_index(vor)]
    #return polys
end

function _voronoi_variance(bon_coordinates, domain_coordinates)
    # Accumulate inclusion weight in each cell, so total = n.
    n = size(bon_coordinates, 2) 
    N = size(domain_coordinates, 2)
    equal_inclusion = n / N

    # Each candidate is assigned to its nearest selected site (Voronoi cell).
    site_tree = KDTree(bon_coordinates)
    cell_totals = zeros(n)
    for i in 1:N
        idxs, _ = knn(site_tree, domain_coordinates[:, i], 1)
        cell_totals[idxs[1]] += equal_inclusion
    end
    return sum((cell_totals .- 1.0) .^ 2) / n
end

function spatialbalance(::VoronoiVariance, domain::Matrix, bon::BiodiversityObservationNetwork)
    bon_coordinates = _FLOAT_TYPE.(bon.coordinates)
    domain_coordinates = _FLOAT_TYPE.(hcat([[ci[1], ci[2]] for ci in CartesianIndices(domain)]...))
    return _voronoi_variance(bon_coordinates, domain_coordinates)    
end

