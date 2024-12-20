




SDMLayer
Matrix

struct Layer{T}
    layer::T
end

struct Stack{T}
    stack::Vector{Layer{T}}
end



abstract type InclusionType end 
struct CategoricalInclusion <: InclusionType
    layer::Layer
    inclusion_probability::Dict
end
struct ContinuousInclusion <: InclusionType
    layer::Layer
end



struct InclusionProbability{<:InclusionType}

end 
