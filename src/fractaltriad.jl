Base.@kwdef struct FractalTriad{I <: Integer, F <: AbstractFloat} <: SpatialSampler
    numsites::I = 50
    horizontal_padding::F = 0.1
    vetical_padding::F = 0.1
    dims::Tuple{I, I}
end

function _generate!(ft::FractalTriad)
    response = zeros(ft.numsites, 2)

    return response
end
