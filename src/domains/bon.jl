"""
    BiodiversityObservationNetwork
"""
struct BiodiversityObservationNetwork{T}
    nodes::Vector{T}
    auxiliary::Matrix
end 

BiodiversityObservationNetwork(nodes::Vector) = BiodiversityObservationNetwork(nodes, missing)

Base.show(io::IO, bon::BiodiversityObservationNetwork) = print(io, "BiodiversityObservationNetwork with $(size(bon)) nodes")

Base.iterate(bon::BiodiversityObservationNetwork) = iterate(bon.nodes)
Base.iterate(bon::BiodiversityObservationNetwork, i) = iterate(bon.nodes, i)

Base.length(bon::BiodiversityObservationNetwork) = length(bon.nodes)
Base.size(bon::BiodiversityObservationNetwork) = length(bon)

Base.push!(bon::BiodiversityObservationNetwork, node) = push!(bon.nodes, node)
Base.getindex(bon::BiodiversityObservationNetwork, i) = bon.nodes[i]

Base.getindex(rs::RasterStack, bon::BiodiversityObservationNetwork) = hcat([rs[n] for n in bon.nodes]...)
Base.getindex(raster::RasterDomain, bon::BiodiversityObservationNetwork) = [raster[n] for n in bon.nodes]



"""
    getpool
"""
getpool(bon::BiodiversityObservationNetwork) = eachindex(bon.nodes)

"""
    getfeatures
"""
getfeatures(bon::BiodiversityObservationNetwork) = bon.auxiliary

"""
    getcoordinates
"""
function getcoordinates(bon::BiodiversityObservationNetwork)
    return Float32.(hcat([[x[1], x[2]] for x in getpool(bon)]...))
end

function BiodiversityObservationNetwork(nodes, raster::RasterDomain{<:SDMLayer})
    Es, Ns = eastings(raster.data), northings(raster.data)
    aux = Matrix([raster[ci] for ci in nodes]')

    nodes = [(Es[n[2]], Ns[n[1]]) for n in nodes]
    return BiodiversityObservationNetwork(nodes, aux)
end

function BiodiversityObservationNetwork(nodes, raster::RasterDomain{<:Matrix})
    aux = Matrix([raster[ci] for ci in nodes]')
    return BiodiversityObservationNetwork(nodes, aux)
end

function BiodiversityObservationNetwork(nodes, rs::RasterStack{<:Matrix})
    aux = hcat([rs[ci] for ci in nodes]...)
    return BiodiversityObservationNetwork(nodes, aux)
end

function BiodiversityObservationNetwork(nodes, rs::RasterStack{<:SDMLayer})
    aux = hcat([rs[ci] for ci in nodes]...)
    Es, Ns = eastings(rs.rasters[1].data), northings(rs.rasters[1].data)
    nodes = [(Es[n[2]], Ns[n[1]]) for n in nodes]
    return BiodiversityObservationNetwork(nodes, aux)
end

function BiodiversityObservationNetwork(nodes, bon::BiodiversityObservationNetwork)
    #aux = hcat([rs[ci] for ci in nodes]...)
    aux = getfeatures(bon)

    return BiodiversityObservationNetwork(nodes, aux)
end
