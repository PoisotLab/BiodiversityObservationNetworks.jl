Base.@kwdef struct FractalTriad{IT<:Integer,FT<:AbstractFloat} <: SpatialSampler 
    numpoints::IT = 50
    padding::FT = 0.1
end

function _generate!(ft::FractalTriad, sdm::M) where {M<:AbstractMatrix}
    response = zeros(ft.numpoints, 2)    

    

    return response
end


