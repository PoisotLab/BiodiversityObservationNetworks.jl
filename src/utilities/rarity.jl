abstract type RarityMetric end 


function _pca(X)
    pca = MVStats.fit(MVStats.PCA, X')
    return MVStats.transform(pca, X')
end 

function _zscore(X)
    z = StatsBase.fit(StatsBase.ZScoreTransform, X)
    return StatsBase.transform(z, X)
end


struct DistanceToMedian <: RarityMetric end 
function rarity(
    ::DistanceToMedian, 
    layers::Vector{<:SDMLayer};
    pca = false
)
    X = hcat([l.grid[l.indices] for l in layers]...)
    X = pca ? _pca(X)' : X
    X = _zscore(X')
    X̄ = map(StatsBase.median, eachrow(X))
    dist = map(xᵢ -> sqrt(sum((xᵢ .- X̄).^2)), eachcol(X))

    rare = deepcopy(first(layers))
    rare.grid[rare.indices] .= dist

    return rare
end

struct MultivariateEnvironmentalSimilarity <: RarityMetric end 


function _mess_score(xᵢⱼ, fᵢⱼ, mⱼ, Mⱼ)
    fᵢⱼ == 0 && return (xᵢⱼ - mⱼ)/(Mⱼ - mⱼ)
    fᵢⱼ == 1 && return (Mⱼ - xᵢⱼ)/(Mⱼ - mⱼ)
    fᵢⱼ < 0.5 && return 2fᵢⱼ
    fᵢⱼ > 0.5 && return 2(1 - fᵢⱼ)
end 

function rarity(
    ::MultivariateEnvironmentalSimilarity,
    layers::Vector{<:SDMLayer}
)
    X = hcat([l.grid[l.indices] for l in layers]...)
    ecdfs = vec(mapslices(StatsBase.ecdf, X; dims = 1))

    mins, maxs = minimum.(layers), maximum.(layers)
    mess = deepcopy(first(layers))

    for (i, cart_idx) in enumerate(findall(mess.indices))
        xᵢ = X[i,:]
        min_Sᵢ = Inf
        for (j, xᵢⱼ) in enumerate(xᵢ) # iterate over each feature
            mⱼ, Mⱼ = mins[j], maxs[j]
            fᵢⱼ = ecdfs[j](xᵢⱼ)
            min_Sᵢ = min(min_Sᵢ, _mess_score(xᵢⱼ, fᵢⱼ, mⱼ, Mⱼ))
        end
        mess[cart_idx] = min_Sᵢ
    end
    return mess
end


