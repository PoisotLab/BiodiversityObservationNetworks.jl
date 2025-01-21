"""
    SimpleRandom

`SimpleRandom` is a [`BONSampler`](@ref) for sampling
[`BiodiversityObservationNetwork`](@ref)s where each location in the spatial
extent has the same probability of inclusion. 
"""
struct SimpleRandom{I<:Integer} <: BONSampler
    number_of_nodes::I
end
_valid_geometries(::SimpleRandom) = (Polygon, Raster, Vector{Polygon}, RasterStack)


# TODO: this should dispatch on individual states if it gets a vector of
# polygons

function _sample(sampler::SimpleRandom, polygons::Vector{<:Polygon}, bon::BiodiversityObservationNetwork)
    vcat([_sample(sampler, poly, bon) for poly in polygons])
end

function _sample(sampler::SimpleRandom, ::Any, bon::BiodiversityObservationNetwork)
    N = sampler.number_of_nodes
    BiodiversityObservationNetwork(Distributions.sample(bon.nodes, N, replace=false))
end 


function _sample(sampler::SimpleRandom, polygon::Polygon)
    x, y = GI.extent(polygon)
    N = sampler.number_of_nodes

    selected_points = Node[]
    ct = 0
    while ct < N
        candidate = (rand(Uniform(x...)), rand(Uniform(y...)))
         if GeometryOps.contains(polygon, candidate)
            push!(selected_points, Node(candidate))
            ct += 1
         end
    end
    return BiodiversityObservationNetwork(selected_points)
end

# Sample w/o replacement from non-empty CIs, then return the long/lat associated
# with it 
function _sample(sampler::SimpleRandom, raster::Raster)
    cart_idxs = nonempty(raster)
    
    node_cidxs = Distributions.sample(cart_idxs, sampler.number_of_nodes, replace=false)

    E,N = SDT.eastings(raster), SDT.northings(raster)

    BiodiversityObservationNetwork([Node((E[I[2]], N[I[1]])) for I in node_cidxs])
end 

_sample(sampler::SimpleRandom, rasters::RasterStack) = _sample(sampler, first(rasters))

function _sample(sampler::SimpleRandom, domain::Vector{<:Polygon})
    @info "You passed a Vector of Polygons."
    @info "Note by default SimpleRandom applies to each Polygon separately"
    @info "To use SimpleRandom on the whole extent, merge the polygons."

    vcat([_sample(sampler, p) for p in domain])
end 


_sample(::SimpleRandom, ::T) where T = throw(ArgumentError("Can't use SimpleRandom on a $T"))
