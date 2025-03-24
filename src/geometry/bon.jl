const __COORDINATE_TYPES = Union{Tuple{<:Real, <:Real}, CartesianIndex}

"""
    Node

A single sampling location within a [`BiodiversityObservationNetwork`](@ref), represented as a coordinate.

A `Node` extends the GeoInterface PointTrait type.
"""
struct Node{C <: __COORDINATE_TYPES}
    coordinate::C
end
Base.show(io::IO, node::Node) = print(io, "📍 Node at $(node.coordinate)")
Base.getindex(node::Node, i) = getindex(node.coordinate, i)

Node(x::T, y::T) where T<:Number = Node((x,y))

GI.isgeometry(::Node)::Bool = true
GI.geomtrait(::Node)::DataType = PointTrait()
GI.ncoord(::PointTrait, ::Node)::Integer = 1
GI.getcoord(::PointTrait, node::Node, i)::Real = node[i]
GI.coordinates(node::Node) = node.coordinate

is_bonifyable(::T) where T = T <: Union{<:__COORDINATE_TYPES,Vector{<:__COORDINATE_TYPES}}


"""
    BiodiversityObservationNetwork

A set of [`Node`](@ref)s that together create a sampling design for monitoring
biodiversity. 

`BiodiversityObservationNetwork` extends the GeoInterface MultiPointTrait type. 
"""
struct BiodiversityObservationNetwork{N <: Node}
    nodes::Vector{N}
end

Base.length(bon::BiodiversityObservationNetwork) = length(bon.nodes)
Base.size(bon::BiodiversityObservationNetwork) = length(bon)
Base.show(io::IO, bon::BiodiversityObservationNetwork) =
    print(io, "🌐 BiodiversityObservationNetwork with $(length(bon)) nodes")
Base.getindex(bon::BiodiversityObservationNetwork, i::Integer) = bon.nodes[i]
Base.getindex(bon::BiodiversityObservationNetwork, i::Vector{<:Integer}) = bon.nodes[i]

Base.iterate(bon::BiodiversityObservationNetwork, i) = iterate(bon.nodes, i)
Base.iterate(bon::BiodiversityObservationNetwork) = iterate(bon.nodes)
Base.vcat(bons::Vector{<:BiodiversityObservationNetwork}) = BiodiversityObservationNetwork(vcat([b.nodes for b in bons]...))


Base.getindex(layer::SDT.SDMLayer, node::Node) = layer[node.coordinate...]
Base.getindex(layer::SDT.SDMLayer, bon::BiodiversityObservationNetwork) = map(c->layer[c], bon)
Base.getindex(layers::Vector{<:SDT.SDMLayer}, bon::BiodiversityObservationNetwork) = hcat([[l[c...] for l in layers] for c in GI.coordinates(bon)]...)



GI.isgeometry(::BiodiversityObservationNetwork)::Bool = true
GI.geomtrait(::BiodiversityObservationNetwork)::DataType = MultiPointTrait()
GI.ngeom(::MultiPointTrait, bon::BiodiversityObservationNetwork)::Integer = length(bon)
GI.getgeom(::MultiPointTrait, bon::BiodiversityObservationNetwork, i) = bon[i]
GI.coordinates(bon::BiodiversityObservationNetwork) = map(GI.coordinates, bon)

function GI.extent(bon::BiodiversityObservationNetwork)
    coords = [n.coordinate for n in bon]
    (bot, top), (left, right) = extrema(first.(coords)), extrema(last.(coords))
    GI.Extent(X=(left, right), Y=(bot, top))
end

features(layers::Vector{<:SDT.SDMLayer}) = findall(first(layers).indices), hcat([[l[x] for l in layers] for x in eachindex(first(layers))]...)

GI.extent(layer::SDMLayer) = begin
    bbox = SDT.boundingbox(layer)
    GI.Extent(X=(bbox[:left], bbox[:right]), Y=(bbox[:bottom], bbox[:top]))
end