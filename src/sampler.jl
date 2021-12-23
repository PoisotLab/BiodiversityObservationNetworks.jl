function Base.rand(method::S, sdm::M) where {S<:SpatialSampler,M<:AbstractMatrix}
    _generate!(method, sdm)
end


