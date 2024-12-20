const __COORDINATE_TYPES = Union{Tuple{<:Real, <:Real}, CartesianIndex}

# TODO: extent this to be GeoInterface.PointLike compat
struct Node{C <: __COORDINATE_TYPES}
    coordinate::C
end
Base.show(io::IO, node::Node) = print(io, "Node at $(node.coordinate)")
Base.getindex(node::Node, i) = getindex(node.coordinate, i)


is_bonifyable(::T) where T = T <: Union{<:__COORDINATE_TYPES,Vector{<:__COORDINATE_TYPES}}


struct BiodiversityObservationNetwork{N <: Node}
    nodes::Vector{N}
end
Base.show(io::IO, bon::BiodiversityObservationNetwork) =
    print(io, "BiodiversityObservationNetwork with $(length(bon)) nodes")
Base.getindex(bon::BiodiversityObservationNetwork, i::Integer) = bon.nodes[i]
Base.length(bon::BiodiversityObservationNetwork) = length(bon.nodes)
Base.iterate(bon::BiodiversityObservationNetwork, i) = iterate(bon.nodes, i)
Base.iterate(bon::BiodiversityObservationNetwork) = iterate(bon.nodes)

Base.vcat(bons::Vector{<:BiodiversityObservationNetwork}) = BiodiversityObservationNetwork(vcat([b.nodes for b in bons]...))
