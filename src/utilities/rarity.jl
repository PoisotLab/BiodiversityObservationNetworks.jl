"""
    RarityMetric

Abstract type encompassing all methods for computing environmental rarity.
"""
abstract type RarityMetric end 


function _zscore(X)
    z = StatsBase.fit(StatsBase.ZScoreTransform, X)
    return z, StatsBase.transform(z, X)
end

function getfeatures(layers; mask=trues(size(first(layers))))
    hcat([[l[i] for l in layers] for i in findall(mask)]...)
end


"""
    DistanceToMedian 

Rarity score defined as Euclidean distance in feature space to the per-feature
median across the raster stack. Optionally, features can be PCA-transformed
prior to z-scoring.
"""
struct DistanceToMedian <: RarityMetric end 
function evaluate(
    ::DistanceToMedian, 
    layers::Vector{<:Matrix};
    mask = trues(size(first(layers)))
)
    X = _FLOAT_TYPE.(getfeatures(layers; mask))
    z, X = _zscore(X)
    X̄ = map(StatsBase.median, eachrow(X))
    dist = map(xᵢ -> sqrt(sum((xᵢ .- X̄).^2)), eachcol(X))

    rare = deepcopy(_FLOAT_TYPE.(first(layers)))
    rare[findall(mask)] .= dist

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

function evaluate(
    ::MultivariateEnvironmentalSimilarity,
    layers::Vector{<:Matrix};
    mask = trues(size(first(layers)))
)
    X = getfeatures(layers; mask)
    ecdfs = vec(mapslices(StatsBase.ecdf, X; dims = 2))


    mins, maxs = minimum.(eachrow(X)), maximum.(eachrow(X))
    mess = deepcopy(_FLOAT_TYPE.(first((layers))))

    pool = findall(mask)

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


"""
    DistanceToAnalogNode <: RarityMetric

TODO
"""
struct DistanceToAnalogNode <: RarityMetric end

"""
    rarity(::DistanceToAnalogNode, bon, layers; pca=false)

For each cell, compute the distance in z-scored feature space to the nearest
selected BON node in `layers`. Optionally apply a shared PCA transform first.
"""
function evaluate(
    ::DistanceToAnalogNode,
    layers::Vector{<:Matrix},
    bon::BiodiversityObservationNetwork;
    mask = trues(size(first(layers)))
)

    X = _FLOAT_TYPE.(getfeatures(layers; mask))
    Xbon = bon.features

    z, X = _zscore(X)
    Xbon = StatsBase.transform(z, Xbon)


    rar = deepcopy(_FLOAT_TYPE.(first(layers)))
    pool = findall(mask)

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

"""
    WithinRange <: RarityMetric

TODO
"""
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
function evaluate(
    ::WithinRange, 
    layers::Vector{<:Matrix},
    bon::BiodiversityObservationNetwork;
    mask = trues(size(first(layers))) 
)
    Xbon = bon.features
    Xextrema = map(extrema, eachrow(Xbon))

    X = getfeatures(layers; mask)
    pool = findall(mask)

    rar = falses(size(first(layers)))
    
    for (i, idx) in enumerate(pool)
        rar[idx] = _point_within_extremes(X[:,i], Xextrema)
    end
    return rar
end 
