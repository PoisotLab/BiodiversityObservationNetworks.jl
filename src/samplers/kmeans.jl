""" 
    KMeans

`KMeans` is a [`BONSampler`](@ref) that generates
[`BiodiversityObservationNetwork`](@ref)s that aim to be representative of the
distribution of environmental variables across an extent represented as a
[`RasterStack`](@ref).

*Arguments*:
- `k`: the number of sampling sites
"""
struct KMeans{I<:Integer} <: BONSampler
    k::I
end 

_valid_geometries(::KMeans) = (RasterStack)


function _sample(::KMeans, ::T) where T 
    @error "Can't use KMeans on a $T"
end 

function _sample(sampler::KMeans, layers::RasterStack)
    cartesian_idx_pool, X = features(layers)
    clusters = Clustering.kmeans(X, sampler.k)

    min_dists = [Inf for _ in 1:sampler.k]
    min_dist_idxs = [0 for _ in 1:sampler.k]

    group_id = clusters.assignments
    μ = clusters.centers
    for i in eachindex(group_id)
        xᵢ = X[:,i]
        gᵢ = group_id[i]
        μᵢ = μ[:,gᵢ]

        dᵢ = sqrt(sum((xᵢ .- μᵢ).^2))
        if dᵢ < min_dists[gᵢ]
            min_dist_idxs[gᵢ] = i 
            min_dists[gᵢ] = dᵢ
        end
    end 

    Es, Ns = eastings(layers), northings(layers)
    
    coords = [(Es[I[2]], Ns[I[1]]) for I in cartesian_idx_pool[min_dist_idxs]] # Cartesian Indices are always in the order (Lat, Long)
    BiodiversityObservationNetwork(Node.(coords))
end 
