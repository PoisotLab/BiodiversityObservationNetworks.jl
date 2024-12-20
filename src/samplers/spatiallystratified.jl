struct SpatiallyStratified{I<:Integer} <: BONSampler
    number_of_nodes::I
end

function _sample(
    sampler::SpatiallyStratified, 
    raster::Raster;
    kwargs...
)
    datatype(raster) <: Integer || throw(ArgumentError("Raster containing spatial strata must be discrete (integer-valued)"))

end



function _sample(
    sampler::SpatiallyStratified, 
    domain::Vector{<:Polygon};
    kwargs...
)
    @info "By default, the number of points within each polygon stratum is proportional to the stratum's area"

    areas = GO.area.(domain)
    areas ./= sum(areas)

    _sample(sampler, domain, areas; kwargs...)
end

function _assign_fixed_inclusions(num_nodes, weights) 
    pts_per_statum = Int.(floor.(weights .* num_nodes))
    leftover = num_nodes - sum(pts_per_statum)
    added_idx = rand(1:length(pts_per_statum), leftover)
    for i in added_idx
        pts_per_statum[i] += 1
    end
    pts_per_statum
end 

function _sample_inclusions(num_nodes, weights) 
    stratum_id = rand(Categorical(weights), num_nodes) 
    pts_per_statum = zeros(Int, length(weights))
    for i in stratum_id
        pts_per_statum[i] += 1
    end
    pts_per_statum
end 

function _sample(
    sampler::SpatiallyStratified, 
    polygons::Vector{<:Polygon}, 
    weights::Vector{<:Real}; 
    fixed_weights = true,
)
    length(polygons) == length(weights) || throw(ArgumentError("Number of provided Polygons is not the same as number of provided weights"))

    N = sampler.number_of_nodes

    nodes_per_stratum = fixed_weights ? _assign_fixed_inclusions(N, weights) : _sample_inclusions(N, weights)
    
    vcat([_sample(SimpleRandom(nodes_per_stratum[i]), polygons[i]) for i in eachindex(polygons) if nodes_per_stratum[i] > 0])
end



