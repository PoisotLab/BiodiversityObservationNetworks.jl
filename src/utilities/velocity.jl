function spatial_gradient(layer)
    offset = CartesianIndices((-1:1, -1:1))
    Δx, Δy = -(SDT.eastings(layer)[[2,1]]...), -(SDT.northings(layer)[[2,1]]...)
    spatial_grad = deepcopy(layer)

    for x in eachindex(layer)
        @inbounds l = layer.grid[x .+ offset]
        @inbounds inc = layer.indices[x .+ offset]

        l[.!(inc)] .= layer.grid[x]
        a,b,c,d,e,f,g,h,i = [l[j,i] for i in 1:3, j in 1:3]

        ∂x = ((c + 2f + i)-(a + 2d + g)) / 8Δx
        ∂y = ((g + 2h + i)-(a + 2b + c)) / 8Δy
        spatial_grad[x] = sqrt((∂x)^2 + (∂y)^2)
    end
    return spatial_grad
end 

function ols(x, y)
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
        _, β = ols(years, y)
        temporal_grad[x] = β
    end
    return temporal_grad
end 

"""
    this follows Loarie et al. 2009 
"""
function velocity(years, timeseries; threshold=0.95)
    sg = spatial_gradient(timeseries[1])
    tg = temporal_gradient(years, timeseries)


    SDT.rescale!.([sg, tg], 1e-4, 1)

    vel = tg / sg
    
    # this ratio gives a very small number of extremely large values, so we clip it to a high percentile that is provided as a kwarg
    τ = Statistics.quantile(vel.grid[vel.indices], [threshold])[1]
    idx = findall(x-> x > τ, vel.grid)
    vel.grid[idx] .= τ
    
    return vel
end



