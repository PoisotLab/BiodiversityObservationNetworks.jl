@kwdef struct BalancedAcceptance{I <: Integer} <: SpatialSampler
    numpoints::I = 50
    α = 1.0
end


function _generate!(piv::BalancedAcceptance, uncertainty::M) where {M<:AbstractMatrix}
    seed = Int32.([floor(10^7*rand()),floor(10^7*rand())])
    np, α = piv.numpoints, piv.α
    x,y = size(uncertainty)

    stduncert = StatsBase.transform(StatsBase.fit(ZScoreTransform, uncertainty, dims=2), uncertainty)
    reluncert = broadcast(x->exp(α*x)/(1+exp(α*x)), stduncert)
    coords = []

    ptct = 0
    while length(coords) < np
        i,j = haltonvalue(seed[1]+ptct,2), haltonvalue(seed[2]+ptct,3)
        candcoord = CartesianIndex(convert.(Int32, [ceil(x*i), ceil(y*j)])...)
        prob = reluncert[candcoord]
        rand() < prob && push!(coords, candcoord)
        ptct += 1
    end
    
    return coords
end

