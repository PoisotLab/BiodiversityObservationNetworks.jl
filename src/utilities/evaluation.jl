"""
    SamplingMetric

Abstract supertype for metrics that evaluate a [`BiodiversityObservationNetwork`](@ref).

Implement `evaluate(metric::MyMetric, result::BiodiversityObservationNetwork)` to add a
new metric.
"""
abstract type SamplingMetric end


"""
    MoransI <: SamplingMetric

Computes Moran's I on the inclusion indicator variable.

# Description
Calculates spatial autocorrelation of the sample indicator ``\\delta`` (1 if sampled, 0 otherwise).
Negative values indicate spatial inhibition (spread), which is desired for balanced sampling.

This version was proposed by [Tille2018MeaSpa](@cite).
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



# ==========================================================================================
# Jensen-Shannon Distance 
# ==========================================================================================

"""
    JensenShannon

The JensenShannon method for evaluating [`BiodiversityObservationNetwork`](@ref) design 
computes how representative the selected BON sites are of the environmental variables 
across the domain using [Jensen-Shannon divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence),
a method for measuring the distance between two distributions. 

# Fields
- `nbins::Int`: number of bins to use when computing empirical probability mass 

"""
@kwdef struct JensenShannon <: SamplingMetric 
    nbins::Int = 25
end

"""
    _standardize

Standardizes the values of a matrix of predictors across the entire population
`Xfull`, and a set of predictors associated with the sampled sites, `Xsampled`,
by scaling each predictor to [0,1].

`Xsampled` is standardized based on the minimum and maximum values of each
predictor across the population, so both matrices return on the same scale.

*Arguments*:
- `Xfull`: an `n` x `d` matrix, where `n` is the size of the population, and `d`
  is the number of predictors
- `Xsampled`: an `m` x `d` matrix, where `m` < `n` is the size of the sample

"""
function _standardize(Xfull, Xsampled)
    Xsampled_std = zeros(size(Xsampled))
    Xfull_std = zeros(size(Xfull))
    for i in axes(Xfull, 1)
        mi, mx = extrema(Xfull[i,:])
        Xfull_std[i,:] .= (Xfull[i,:] .- mi) ./ (mx - mi)
        Xsampled_std[i,:] = (Xsampled[i,:] .- mi) ./ (mx - mi)
    end 
    return Xfull_std, Xsampled_std
end

"""
    _histbin(x, edges)

Compute empirical probability mass by binning empirical values for a single feature `x`.
"""
function _histbin(x, edges)
    n = length(edges) - 1
    counts = fill(1e-10, n)  # avoids log(0) in KL divergence
    for v in x
        b = clamp(searchsortedlast(edges, v), 1, n)
        counts[b] += 1
    end
    return counts ./ sum(counts)
end

"""
    _jensen_shannon(Xfull, Xsampled; nbins=20)

Compute Jensen-Shannon divergence between two empirical distributions of different size:
`Xfull`, the features across the entire domain, and `Xsampled`, the features at sampled sites.
"""
function _jensen_shannon(Xfull, Xsampled; nbins=20)
    Xf, Xs = _standardize(Xfull, Xsampled)
    bin_edges = range(0.0, 1.0, length=nbins + 1)
    d = size(Xf, 1)
    total = 0.0

    for i in 1:d
        # P is the full distribution of the feature across domain
        P = _histbin(view(Xf, i, :), bin_edges)

        # Q is distribution of feature at sampled sites
        Q = _histbin(view(Xs, i, :), bin_edges)

        # Mixture distribution of P and Q
        M = 0.5 .* (P .+ Q)
        
        # Definition of Jensen-Shannon from KL-Divergence 
        total += 0.5 * StatsBase.kldivergence(P, M) + 0.5 * StatsBase.kldivergence(Q, M)
    end

    # return average JS across each dimension
    return total / d
end

""" 
    evaluate(::JensenShannon, bon::BiodiversityObservationNetwork, layers)
"""
function evaluate(
    js::JensenShannon,
    layers::Vector{<:Matrix}, 
    bon::BiodiversityObservationNetwork,
)

    Xfull = hcat([[l[i] for l in layers] for i in findall(x -> x isa Number, layers[begin])]...)
    Xbon = hcat([[l[i] for l in layers] for i in bon.sites]...)

    return _jensen_shannon(Xfull, Xbon; nbins = js.nbins)
end 