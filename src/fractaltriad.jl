@kwdef struct FractalTriad{IT<:Integer,FT<:AbstractFloat} <: SpatialSampler 
    numpoints::IT = 50
    padding::FT = 0.1
end

function _generate!(response::M, ft::FractalTriad, sdm::M) where {M<:AbstractMatrix}
    @assert size(response) == size(sdm)


end


