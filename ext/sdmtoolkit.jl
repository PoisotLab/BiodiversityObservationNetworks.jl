module SDTExt

using BiodiversityObservationNetworks
using SpeciesDistributionToolkit

@info "Loading BONs.jl support for SimpleSDMLayers.jl ..."

function stack(
    layers::Vector{<:SpeciesDistributionToolkit.SimpleSDMLayers.SimpleSDMLayer},
)
    # assert all layers are the same size if we add this function to BONs.jl
    mat = zeros(size(first(layers))..., length(layers))
    for (l, layer) in enumerate(layers)
        thismin, thismax = extrema(layer)
        mat[:, :, l] .= broadcast(
            x -> isnothing(x) ? NaN : (x - thismin) / (thismax - thismin),
            layer.grid,
        )
    end
    return mat
end

end
