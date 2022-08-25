"""
    entropize!(U::Matrix{AbstractFloat}, A::Matrix{Number})

This function turns a matrix `A` (storing measurement values) into pixel-wise
entropy values, stored in a matrix `U` (that is previously allocated).

Pixel-wise entropy is determined by measuring the empirical probability of
randomly picking a value in the matrix that is either lower or higher than the
pixel value. The entropy of both these probabilities are calculated using the
-p×log(2,p) formula. The entropy of the pixel is the *sum* of the two entropies,
so that it is close to 1 for values close to the median, and close to 0 for
values close to the extreme of the distribution.
"""
function entropize!(U::Matrix{RT}, A::Matrix{T}) where {RT <: AbstractFloat, T <: Number}
    @assert basis > zero(basis)
    for i in eachindex(A)
        p_high = mean(A .< A[i])
        p_low = mean(A .>= A[i])
        e_high = p_high .* log2(p_high)
        e_low = p_low .* log2(p_low)
        U[i] = -e_high .- e_low
    end
    U[findall(isnan, U)] .= zero(eltype(U))
    m = maximum(U)
    for i in eachindex(U)
            U[i] /= m
        end
    return U
end

"""
    entropize(A::Matrix{Number})

Allocation version of `entropize!`.
"""
function entropize(A::Matrix{T}) where {T <: Number}
    U = zeros(size(A))
    return entropize!(U, A)
end