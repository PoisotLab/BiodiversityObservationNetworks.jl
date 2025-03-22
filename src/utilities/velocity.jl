abstract type VelocityMetric end 

struct Loarie2009 <: VelocityMetric end 
struct ClosestAnalogue <: VelocityMetric end 

function spatial_gradient(layer)
    offset = CartesianIndices((-1:1, -1:1))
    Δx, Δy = -(SDT.eastings(layer)[[2,1]]...), -(SDT.northings(layer)[[2,1]]...)
    spatial_grad = deepcopy(layer.raster)

    for x in eachindex(layer)
        @inbounds l = layer.raster.grid[x .+ offset]
        @inbounds inc = layer.raster.indices[x .+ offset]

        l[.!(inc)] .= layer.raster.grid[x]
        a,b,c,d,e,f,g,h,i = [l[j,i] for i in 1:3, j in 1:3]

        ∂x = ((c + 2f + i)-(a + 2d + g)) / 8Δx
        ∂y = ((g + 2h + i)-(a + 2b + c)) / 8Δy
        spatial_grad[x] = sqrt((∂x)^2 + (∂y)^2)
    end
    return spatial_grad
end 

function _ols(x, y)
    X = hcat(ones(size(x, 1)), x) 
    XᵀX⁻¹ = inv(X' * X)
    α, β = XᵀX⁻¹ * X' * y
    return α, β
end

function temporal_gradient(years, timeseries)
    baseline = first(timeseries)
    temporal_grad = deepcopy(baseline)
    for x in eachindex(baseline)
        y = [l[x] for l in timeseries]
        _, β = _ols(years, y)
        temporal_grad[x] = β
    end
    return temporal_grad
end 

"""
    this follows Loarie et al. 2009 
"""
velocity(::Type{Loarie2009}, args...) = velocity(Loarie2009(), args...)
function velocity(::Loarie2009, years, timeseries; threshold=0.95)
    sg = spatial_gradient(timeseries[1])
    tg = temporal_gradient(years, timeseries)
    vel = tg / sg
    
    # this ratio gives a very small number of extremely large values, so we clip it to a high percentile that is provided as a kwarg
    τ = Statistics.quantile(vel.grid[vel.indices], [threshold])[1]
    idx = findall(x-> x > τ, vel.grid)
    vel.grid[idx] .= τ
    
    return vel
end



function _nearest_feature_neighbor(baseline, future)
    cart_idx, baseline_features = features(baseline)
    _, future_features = features(future)

    # features should be zscored so units are in SD, otherwise distance is relative to units for each feature

    closest_analogue = fill(CartesianIndex(0,0), size(baseline))

    kd = NearestNeighbors.KDTree(future_features)
    for (i,bi) in enumerate(eachcol(baseline_features))
        nearest_idx, _ = NearestNeighbors.knn(kd, bi, 1)
        closest_analogue[cart_idx[i]] = cart_idx[nearest_idx[begin]]
    end
    return closest_analogue
end

_euclidian_dist(x,y) = sqrt(sum((x .- y).^2))

velocity(::Type{ClosestAnalogue}, args...) = velocity(ClosestAnalogue(), args...)
function velocity(::ClosestAnalogue, baseline, future)
    closest_analogues = _nearest_feature_neighbor(baseline, future)

    Es, Ns = SDT.eastings(baseline), SDT.northings(baseline)

    vel = baseline isa RasterStack ? deepcopy(baseline[1].raster) : deepcopy(baseline.raster)
    for ci in findall(vel.indices)
        future_cart = closest_analogues[ci]
        base_coord = Es[ci[2]], Ns[ci[1]]
        future_coord = Es[future_cart[2]], Ns[future_cart[1]]
        vel[ci] = _euclidian_dist(base_coord, future_coord)
    end
    return vel
end 
