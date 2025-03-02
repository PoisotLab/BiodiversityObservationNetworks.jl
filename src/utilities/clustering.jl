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

abstract type ClusteringAlgorithm end 
@kwdef struct KMeans{I<:Integer} <: ClusteringAlgorithm 
    k::I = 3
end

@kwdef struct KMedoids{I<:Integer} <: ClusteringAlgorithm 
    k::I = 3
end

function cluster!(km::KMeans, layers::RasterStack, assignments::Raster)
    k = km.k 
    idx, mat = features(layers)
    clust = Clustering.kmeans(mat, k)
    assignments.raster.grid[idx] .= clust.assignments
    return assignments
end 


function cluster!(km::KMedoids, layers::RasterStack, assignments::Raster)
    k = km.k 
    idx, mat = features(layers)
    dist_mat = _feature_distance(mat)
    clust = Clustering.kmedoids(dist_mat, k)
    assignments.raster.grid[idx] .= clust.assignments
    return assignments
end 

function cluster(alg::ClusteringAlgorithm, layers::RasterStack)
    assignments = deepcopy(first(layers))
    cluster!(alg, layers, assignments)
end