"""
abstract type ClusteringAlgorithm end
function cluster(alg::ClusteringAlgorithm, layers::Vector{<:SDMLayer})
    assignments = deepcopy(first(layers))
    cluster!(alg, layers, assignments)
end


@kwdef struct KMeans{I<:Integer} <: ClusteringAlgorithm 
    k::I = 3
end

function cluster(km::KMeans, stack::Vector{<:SDMLayer}; kw...)
    k = km.k 
    layers = map(x->x.raster, stack.stack)
    clustering_result = Clustering.kmeans(layers, k; kw...)
    return SDMLayer(clustering_result, layers)
end 


@kwdef struct FuzzyCMeans{I<:Integer,F<:Real} <: ClusteringAlgorithm 
    c::I = 3
    fuzzyness::F = 2.0
end

function cluster(cm::FuzzyCMeans, stack::Vector{<:SDMLayer}; kw...)
    c, fuzzyness = cm.c, cm.fuzzyness
    layers = map(x->x.raster, stack.stack)
    clustering_result = Clustering.fuzzy_cmeans(layers, c, fuzzyness; kw...)
    return SDMLayer(clustering_result, layers)
end 

"""
