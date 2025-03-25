abstract type RarityMetric end 

function _pca(X)
    pca = MVStats.fit(MVStats.PCA, X')
    return MVStats.transform(pca, X')
end 

function _zscore(X)
    z = StatsBase.fit(StatsBase.ZScoreTransform, X)
    return z, StatsBase.transform(z, X)
end


struct DistanceToMedian <: RarityMetric end 
function rarity(
    ::DistanceToMedian, 
    layers::Vector{<:SDMLayer};
    pca = false
)
    X = hcat([l.grid[l.indices] for l in layers]...)
    X = pca ? _pca(X)' : X
    z, X = _zscore(X')
    X̄ = map(StatsBase.median, eachrow(X))
    dist = map(xᵢ -> sqrt(sum((xᵢ .- X̄).^2)), eachcol(X))

    rare = deepcopy(first(layers))
    rare.grid[rare.indices] .= dist

    return rare
end

"""
    MultivariateEnvironmentalSimilarity

Multivariate-Environmental Similarity Score (MESS) is a metric introduced by [Elith2010ArtModelling](@cite) to quantify environmental similarity, meaning _lower_ values indicate _more rare_ environments. It is conceptually similar to how BioClim scores work, but enables negative values. 
"""
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


struct DistanceToAnalogNode <: RarityMetric end
function rarity(
    ::DistanceToAnalogNode, 
    layers::Vector{<:SDMLayer}, 
    bon::BiodiversityObservationNetwork;
    pca = false
)
    X = hcat([l.grid[l.indices] for l in layers]...)
    X = pca ? _pca(X)' : X
    z, X = _zscore(X')

    Xbon = StatsBase.transform(z, layers[bon])

    rarity = deepcopy(first(layers))

    for (i, ci) in enumerate(eachindex(first(layers)))
        Xi = X[:,i]
        min_dist = Inf
        for Xb in eachcol(Xbon)
            dist = sqrt(sum((Xb .- Xi).^2))
            min_dist = dist < min_dist ? dist : min_dist
        end 
        rarity[ci] = min_dist
    end
    return rarity
end 
