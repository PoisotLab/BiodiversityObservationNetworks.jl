@kwdef struct AdaptiveSpatialSampling{I <: Integer} <: SpatialSampler
    numpoints::I = 50
end

function _generate!(ass::AdaptiveSpatialSampling, uncertainty::M) where {M<:AbstractMatrix}
    
    return coords
end
