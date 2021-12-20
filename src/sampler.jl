import Random.rand!

function Base.rand(method::S, sdm::M) where {S<:SpatialSampler,M<:AbstractMatrix}
    
end

function rand!(mat::M where N, method::T, sdm::M) where {T<:SpatialSampler,M<:AbstractMatrix}
    _generate!(mat, method, sdm)
end

