function stack(layers::Vector{<:AbstractMatrix})
    # assert all layers are the same size if we add this function to BONs.jl

    mat = zeros(size(first(layers))..., length(layers))
    for (l, layer) in enumerate(layers)
        nonnanvals = vec(layer[findall(!isnan, layer)])
        thismin, thismax = findmin(nonnanvals)[1], findmax(nonnanvals)[1]

        mat[:, :, l] .=
            broadcast(x -> isnan(x) ? NaN : (x - thismin) / (thismax - thismin), layer)
    end
    return mat
end
