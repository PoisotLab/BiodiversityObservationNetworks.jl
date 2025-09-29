"""
    RarityMetric

Abstract type encompassing all methods for computing environmental rarity.
"""
abstract type RarityMetric end 


rarity(rm::RarityMetric, bon::BiodiversityObservationNetwork, domain; kwargs...) = rarity(rm, bon, to_domain(domain); kwargs...)

rarity(rm::RarityMetric, domain; kwargs...) = rarity(rm, to_domain(domain); kwargs...)


function _fit_pca(X)
    return MultivariateStats.fit(MultivariateStats.PCA, X)
end 
function _transform_pca(pca, X)
    return MultivariateStats.transform(pca, X)
end 

function _zscore(X)
    z = StatsBase.fit(StatsBase.ZScoreTransform, X)
    return z, StatsBase.transform(z, X)
end


"""
    DistanceToMedian 

Rarity score defined as Euclidean distance in feature space to the per-feature
median across the raster stack. Optionally, features can be PCA-transformed
prior to z-scoring.
"""
struct DistanceToMedian <: RarityMetric end 
function rarity(
    ::DistanceToMedian, 
    layers::RasterStack;
    pca = false
)
    X = getfeatures(layers)
    if pca 
        X = _transform_pca(_fit_pca(X), X)
    end
    z, X = _zscore(X)
    X̄ = map(StatsBase.median, eachrow(X))
    dist = map(xᵢ -> sqrt(sum((xᵢ .- X̄).^2)), eachcol(X))


    pool = getpool(layers)

    rare = deepcopy(first(layers))
    rare.data.grid[pool] .= dist

    return rare
end

"""
    MultivariateEnvironmentalSimilarity

Multivariate Environmental Similarity Surface (MESS). For each cell, compute the
minimum over features of a per-feature similarity score derived from the ECDF of
the training distribution, following the standard MESS definition.
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
    layers::RasterStack
)
    X = getfeatures(layers)
    ecdfs = vec(mapslices(StatsBase.ecdf, X; dims = 2))

    mins, maxs = minimum.(layers.rasters), maximum.(layers.rasters)
    mess = deepcopy(first(layers))

    pool = getpool(layers)

    for (i, cart_idx) in enumerate(pool)
        xᵢ = X[:,i]
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

"""
    rarity(::DistanceToAnalogNode, bon, layers; pca=false)

For each cell, compute the distance in z-scored feature space to the nearest
selected BON node in `layers`. Optionally apply a shared PCA transform first.
"""
function rarity(
    ::DistanceToAnalogNode, 
    bon::BiodiversityObservationNetwork,
    layers::RasterStack;
    pca = false
)

    X = getfeatures(layers)
    Xbon = layers[bon]
    if pca
        pca_fit = _fit_pca(X)
        X = _transform_pca(pca_fit, X)
        Xbon = _transform_pca(pca_fit, Xbon)
    end

    z, X = _zscore(X)
    Xbon = StatsBase.transform(z, Xbon)


    rar = deepcopy(first(layers))
    pool = getpool(layers)

    for (i, ci) in enumerate(pool)
        Xi = X[:,i]
        min_dist = Inf
        for Xb in eachcol(Xbon)
            dist = sqrt(sum((Xb .- Xi).^2))
            min_dist = dist < min_dist ? dist : min_dist
        end 
        rar[ci] = min_dist
    end
    return rar
end 


struct WithinRange <: RarityMetric end

function _point_within_extremes(point, extremes)
    for (j, xᵢ) in enumerate(point)
        extremes[j][1] <= xᵢ <= extremes[j][2] || return false
    end
    return true
end     

"""
    rarity(::WithinRange, bon, layers)

Boolean rarity surface indicating whether each cell lies within the hyper-
rectangle spanned by the per-feature minima and maxima of the BON nodes.
"""
function rarity(
    ::WithinRange, 
    bon::BiodiversityObservationNetwork,
    layers::RasterStack, 
)
    Xbon = layers[bon]
    Xextrema = map(extrema, eachrow(Xbon))

    X = getfeatures(layers)
    pool = getpool(layers)

    rar = deepcopy(first(layers))
    
    for (i, idx) in enumerate(pool)
        rar[idx] = _point_within_extremes(X[:,i], Xextrema)
    end
    return rar
end 
