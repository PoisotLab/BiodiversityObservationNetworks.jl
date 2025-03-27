"""
    SpatiallyStratified

`SpatiallyStratified` is a [`BONSampler`](@ref) for choosing sites across a set of different spatial stratum. 
"""
@kwdef struct SpatiallyStratified{I<:Integer} <: BONSampler
    number_of_nodes::I = 30
end

function _sample(
    sampler::SpatiallyStratified, 
    raster::SDMLayer{T},
    args...;
    kwargs...
) where T
    T <: Integer || throw(ArgumentError("Raster containing spatial strata must be discrete (integer-valued)"))
end
 
function _sample(
    sampler::SpatiallyStratified,
    raster::SDMLayer;
    kwargs...
)
    @info "By default, the number of points within each stratum is proportional to the stratum's area"

    strata = unique(raster)
    areas = length.([findall(isequal(s), raster) for s in strata]) ./ sum(raster.indices)
    _sample(sampler, raster, areas; kwargs...)
end

function _sample(
    sampler::SpatiallyStratified,
    raster::SDMLayer,
    weights::Vector{<:Real}; 
    fixed_weights = true,
)
    strata = unique(raster)
    idx_per_strata = [findall(isequal(s), raster) for s in strata]
    nodes_per_stratum = _get_nodes_per_stratum(sampler, weights, fixed_weights)

    selected_idx = vcat([Distributions.sample(idx_per_strata[i], n, replace=false) for (i,n) in enumerate(nodes_per_stratum)]...)

    Es, Ns = SDT.eastings(raster), SDT.northings(raster)

    return BiodiversityObservationNetwork([Node(Es[i[2]], Ns[i[1]]) for i in selected_idx])
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

function _get_nodes_per_stratum(sampler, weights, fixed_weights)
    N = sampler.number_of_nodes
    fixed_weights ? _assign_fixed_inclusions(N, weights) : _sample_inclusions(N, weights)
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


    nodes_per_stratum = _get_nodes_per_stratum(sampler,weights, fixed_weights)
    
    

    # todo: there's no reason this can't be a sampler from a set of acceptable
    # samplers
    # this proves harder when it doesn't dispatch on the same tpe that
    # SpatialStratified does, but it may not be impossible. 
        # if it has the same dispatch its ez. 

    vcat([_sample(SimpleRandom(nodes_per_stratum[i]), polygons[i]) for i in eachindex(polygons) if nodes_per_stratum[i] > 0])
end

# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use SpatiallyStratified with default constructor on a vector of Polygons" begin
    #=polys = openstreetmap.(["France", "Germany", "Belgium"])
    ss = SpatiallyStratified()
    bon = sample(ss, polys)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == ss.number_of_nodes=#
end


@testitem "We can use SpatiallyStratified with default constructor on a discrete Raster" begin
    mat = zeros(Int, 10,10)
    mat[6:end, 5:end] .= 1
    mat[1:5,5:end] .= 2

    raster = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(mat)
    ss = SpatiallyStratified()
    bon = sample(ss, raster)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == ss.number_of_nodes
end
