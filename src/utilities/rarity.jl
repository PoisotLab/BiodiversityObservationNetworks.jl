abstract type RarityMetric end 

struct DistanceToMedian <: RarityMetric end 

function _pca(X)
    pca = MVStats.fit(MVStats.PCA, X')
    return MVStats.transform(pca, X')
end 

function rarity(
    ::DistanceToMedian, 
    layers::Vector{<:SDMLayer};
    pca = false
)
    X = hcat([l.grid[l.indices] for l in layers]...)
    X = pca ? _pca(X)' : X
    z = StatsBase.fit(StatsBase.ZScoreTransform, X')
    zfeat = StatsBase.transform(z, X')
    X̄ = map(StatsBase.median, eachrow(zfeat))
    dist = map(xᵢ -> sqrt(sum((xᵢ .- X̄).^2)), eachcol(zfeat))

    rare = deepcopy(first(layers))
    rare.grid[rare.indices] .= dist

    return rare
end


