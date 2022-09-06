@info "Loading BONs.jl support for SimpleSDMLayers.jl ..."


function BiodiversityObservationNetworks.stack(layers::Vector{<:SimpleSDMLayers.SimpleSDMLayer})
    # assert all layers are the same size if we add this function to BONs.jl
    mat = zeros(size(first(layers))..., length(layers))
    for (l,layer) in enumerate(layers)
        mat[:,:,l] .= broadcast(x->isnothing(x) ? NaN : x, layer.grid)
    end
    mat
end
