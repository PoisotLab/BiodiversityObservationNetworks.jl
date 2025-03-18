#= 
function _feature_distance(feat)
    nfeat = size(feat, 2)
    feat_dist = zeros(nfeat, nfeat)
    for i in 1:nfeat, j in i+1:nfeat
        dist = sqrt(sum((feat[:,i] .- feat[:,j]).^2)) 
        feat_dist[i,j] = dist
        feat_dist[j,i] = dist
    end 
    return feat_dist
end
=#

abstract type ClusteringAlgorithm end
function cluster(alg::ClusteringAlgorithm, layers::RasterStack)
    assignments = deepcopy(first(layers))
    cluster!(alg, layers, assignments)
end


@kwdef struct KMeans{I<:Integer} <: ClusteringAlgorithm 
    k::I = 3
end

function cluster(km::KMeans, stack::RasterStack; kw...)
    k = km.k 
    layers = map(x->x.raster, stack.stack)
    clustering_result = Clustering.kmeans(layers, k; kw...)
    return Raster(SDMLayer(clustering_result, layers))
end 


@kwdef struct FuzzyCMeans{I<:Integer,F<:Real} <: ClusteringAlgorithm 
    c::I = 3
    fuzzyness::F = 2.0
end

function cluster(cm::FuzzyCMeans, stack::RasterStack; kw...)
    c, fuzzyness = cm.c, cm.fuzzyness
    layers = map(x->x.raster, stack.stack)
    clustering_result = Clustering.fuzzy_cmeans(layers, c, fuzzyness; kw...)
    return RasterStack(SDMLayer(clustering_result, layers))
end 


