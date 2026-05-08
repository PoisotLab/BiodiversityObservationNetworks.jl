"""
    VelocityMetric

Abstract type encompassing all methods for computing environmental velocity.
"""
abstract type VelocityMetric end 

# ---------------------------------------------
# Begin Loarie2009 

struct Loarie2009 <: VelocityMetric end 

spatial_gradient(layer) = nothing

function _ols(x, y)
    X = hcat(ones(size(x, 1)), x) 
    XᵀX⁻¹ = inv(X' * X)
    α, β = XᵀX⁻¹ * X' * y
    return α, β
end

function temporal_gradient(years, timeseries)
    baseline = first(timeseries)
    temporal_grad = _FLOAT_TYPE.(copy(baseline))
    for x in eachindex(baseline)
        y = [l[x] for l in timeseries]
        _, β = _ols(years, y)
        temporal_grad[x] = β
    end
    return temporal_grad
end 


velocity(::Type{Loarie2009}, args...) = velocity(Loarie2009(), args...)

# End Loarie2009 
# ---------------------------------------------


# ---------------------------------------------
# Begin ClosestAnalogue

struct ClosestAnalogue <: VelocityMetric end 

function _nearest_feature_neighbor(baseline, future)
    cart_idx, baseline_features = getfeatures(baseline)
    _, future_features = getfeatures(future)

    # features should be zscored so units are in SD, otherwise distance is relative to units for each feature
    closest_analogue = fill(CartesianIndex(0,0), size(first(baseline)))

    kd = NearestNeighbors.KDTree(future_features)
    for (i,bi) in enumerate(eachcol(baseline_features))
        nearest_idx, _ = NearestNeighbors.knn(kd, b, 1)
        closest_analogue[cart_idx[i]] = cart_idx[nearest_idx[begin]]
    end
    return closest_analogue
end

_euclidian_dist(x,y) = sqrt(sum((x .- y).^2))

velocity(::Type{ClosestAnalogue}, args...) = velocity(ClosestAnalogue(), args...)
#velocity(::ClosestAnalogue, baseline::SDMLayer, future::SDMLayer) = velocity(ClosestAnalogue(), [baseline], [future])
 
function velocity(::ClosestAnalogue, baseline, future)
    closest_analogues = _nearest_feature_neighbor(baseline, future)
    Es = SDT.eastings(first(baseline))
    Ns = SDT.northings(first(baseline))
    vel = deepcopy(baseline[1]) 
    for ci in findall(vel.indices)
        future_cart = closest_analogues[ci]
        base_coord = Es[ci[2]], Ns[ci[1]]
        future_coord = Es[future_cart[2]], Ns[future_cart[1]]
        vel[ci] = _euclidian_dist(base_coord, future_coord)
    end
    return vel
end 

# End ClosestAnalogue 
# ---------------------------------------------
