"""
    SimpleRandom

`SimpleRandom` is a [`BONSampler`](@ref) for sampling
[`BiodiversityObservationNetwork`](@ref)s where each location in the spatial
extent has the same probability of inclusion. 
"""
@kwdef struct SimpleRandom{I<:Integer} <: BONSampler
    number_of_nodes::I = _DEFAULT_NUM_NODES
end
_valid_geometries(::SimpleRandom) = (Polygon, Raster, Vector{Polygon}, RasterStack, BiodiversityObservationNetwork)


# TODO: this should dispatch on individual states if it gets a vector of
# polygons
function _sample(sampler::SimpleRandom, polygons::Vector{<:Polygon}, bon::BiodiversityObservationNetwork)
    vcat([_sample(sampler, poly, bon) for poly in polygons])
end

# For multistage
function _sample(sampler::SimpleRandom, ::Any, bon::BiodiversityObservationNetwork)
    N = sampler.number_of_nodes
    BiodiversityObservationNetwork(Distributions.sample(bon.nodes, N, replace=false))
end 


function _sample(sampler::SimpleRandom, bon::BiodiversityObservationNetwork)
    N = sampler.number_of_nodes
    BiodiversityObservationNetwork(Distributions.sample(bon.nodes, N, replace=false))
end 


function _sample(sampler::SimpleRandom, polygon::Polygon)
    x, y = GI.extent(polygon)
    _londist, _latdist = Uniform(x...), Uniform(y...)
    N = sampler.number_of_nodes

    selected_points = Node[]
    ct = 0
    while ct < N
        candidate = (rand(_londist), rand(_latdist))
         if GeometryOps.contains(polygon, candidate)
            push!(selected_points, Node(candidate))
            ct += 1
         end
    end
    return BiodiversityObservationNetwork(selected_points)
end

# Sample w/o replacement from non-empty CIs, then return the long/lat associated
# with it 
function _sample(sampler::SimpleRandom, raster::SDMLayer)
    cart_idxs = findall(raster.indices)
    
    node_cidxs = Distributions.sample(cart_idxs, sampler.number_of_nodes, replace=false)

    E,N = SDT.eastings(raster), SDT.northings(raster)

    BiodiversityObservationNetwork([Node((E[I[2]], N[I[1]])) for I in node_cidxs])
end 

_sample(sampler::SimpleRandom, rasters::Vector{<:SDMLayer}) = _sample(sampler, first(rasters))

function _sample(sampler::SimpleRandom, domain::Vector{<:Polygon})
    @info "You passed a Vector of Polygons."
    @info "Note by default SimpleRandom applies to each Polygon separately"
    @info "To use SimpleRandom on the whole extent, merge the polygons."

    vcat([_sample(sampler, p) for p in domain])
end 


_sample(::SimpleRandom, ::T) where T = throw(ArgumentError("Can't use SimpleRandom on a $T"))

# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use SimpleRandom with default arguments on a Raster" begin
    raster = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(zeros(50, 100))
    srs = SimpleRandom()
    bon = sample(srs, raster)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == srs.number_of_nodes
end


@testitem "We can use SimpleRandom with default arguments on a Polygon" begin
    poly = openstreetmap("COL")
    srs = SimpleRandom()
    bon = sample(srs, poly)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == srs.number_of_nodes
end