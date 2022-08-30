"""
    BalancedAcceptance

A `BONRefiner`.



"""
Base.@kwdef mutable struct Uniqueness{I <: Integer} <: BONRefiner
    numpoints::I = 30
    function Uniqueness(numpoints)
        if numpoints < one(numpoints)
            throw(
                ArgumentError(
                    "You cannot have a Uniqueness sampler with less than one point",
                ),
            )
        end
        return new{typeof(numpoints)}(numpoints)
    end
end


function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::Uniqueness,
    layers::Array{T,N},
) where {T <: AbstractFloat,N}
    
    N <= 2 && throw(ArgumentError("Uniqueness needs more than one layer to work."))

    covscore = zeros(length(pool))
    for (i,p1) in enumerate(pool)
        v1 = layers[p1[1], p1[2], :]
        for (j,p2) in enumerate(pool)
            v2 = layers[p2[1], p2[2],:]
            if p1 != p2
                covscore[i] += abs(cov(v1,v2))
            end
        end 
    end

    np = sampler.numpoints
    sortedvals = sortperm(vec(covscore))
    
    coords[:] .= pool[sortedvals[1:np]]
    return (coords, layers)
end