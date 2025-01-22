# NOTE: Order of auxilliary variables is important for cube sampling because 
# if a balanced sample is not possible, they are relaxed in succession from last
# variable
# would we ever want to use the inefficient method? 
    # its O(n²) as opposed to O(n), where n is population size 

"""
    CubeSampling

"""
Base.@kwdef struct CubeSampling{I<:Integer} <: BONSampler
    number_of_nodes::I = 50
end 

_valid_geometries(::CubeSampling) = (BiodiversityObservationNetwork, RasterStack, Raster)


# for multistage
function _sample(sampler::CubeSampling, layers::RasterStack, bon::BiodiversityObservationNetwork)
    cart_idx = _get_cartesian_idx(layers, bon)
    feat = layers[cart_idx]
    N = sampler.number_of_nodes
    π_optimal, candidate_pool = _cube(N, cart_idx, feat)
    Es, Ns = eastings(layers), northings(layers)
    selected = candidate_pool[findall(isequal(1), π_optimal)]
    return BiodiversityObservationNetwork([Node(Es[idx[2]], Ns[idx[1]]) for idx in selected])

end

function _sample(sampler::CubeSampling, layers::RasterStack)
    cart_idx, feat = features(layers)
    N = sampler.number_of_nodes
    π_optimal, candidate_pool = _cube(N, cart_idx, feat)
    Es, Ns = eastings(layers), northings(layers)
    selected = candidate_pool[findall(isequal(1), π_optimal)]
    return BiodiversityObservationNetwork([Node(Es[idx[2]], Ns[idx[1]]) for idx in selected])
end

function _cube(N, cart_idx, feat)
    πₖ = [N / size(feat, 2) for _ in axes(feat, 2)]

    # sort points by distance in auxillary variable space
    dist = mahalanobis(πₖ, feat)
    
    # high to low
    perm = sortperm(dist) |> reverse
    
    candidate_pool, πₖ = cart_idx[perm], πₖ[perm]

    x = feat[:, perm]
    x = vcat(transpose(πₖ), x) # add inclusion probabilities to fix sample size

    π_optimal_flight = cubefastflight(πₖ, x)
    # check if there are non-integer probabilities
    non_int_idx = findall(u -> u != 0 && u != 1, π_optimal_flight)
    # if so, perform landing phase to resolve them
    π_optimal = isempty(non_int_idx) ? π_optimal_flight : cubeland(π_optimal_flight, πₖ, x)

    return π_optimal, candidate_pool
end


#TODO: refactor
function mahalanobis(πₖ, x)
    # drop variables that are the same for all points
    num_uniq = map(y -> length(unique(y)), eachrow(x))
    nonuniq_ind = findall(z -> z == 1, num_uniq)
    x = length(nonuniq_ind) > 0 ? x[1:end .!= nonuniq_ind, :] : x

    N = length(πₖ)
    p = size(x, 1)

    x̂ = x ./ reshape(πₖ, :, N)
    mean_x̂ = (1 / N) .* sum(x̂; dims = 2)

    k_vecs = x̂ .- mean_x̂
    outer_prods = Array{Float64}(undef, p, p, N)
    for i in 1:N
        outer_prods[:, :, i] = k_vecs[1:p, i] * transpose(k_vecs[1:p, i])
    end

    sigma = (1 / (N - 1)) * dropdims(sum(outer_prods; dims = 3); dims = 3)
    inv_sigma = inv(sigma)

    d = Vector{Float64}(undef, N)
    for i in 1:N
        d[i] = (transpose(x̂[1:p, i] - mean_x̂) * inv_sigma * (x̂[1:p, i] - mean_x̂))[1]
    end

    return d
end 

function cubefastflight(πₖ, x)

    # number of auxillary variables
    num_aux = size(x)[1]

    # get all non-integer probabilities
    non_int_ind = findall(u -> u != 0 && u != 1, πₖ)

    π_nonint = πₖ[non_int_ind]

    Ψ = π_nonint[1:(num_aux + 1)]
    r = collect(1:(num_aux + 1))

    A = x[:, non_int_ind] ./ reshape(π_nonint, :, length(non_int_ind))
    B = A[:, 1:(num_aux + 1)]

    k = num_aux + 2

    while k <= length(π_nonint)
        Ψ = update_psi(Ψ, B)

        if length(findall(z -> z .< 0, Ψ)) > 0
            throw(error("Negative inclusion probability"))
        end

        # update for the probabilities that are now integers
        i = 0
        while i < length(Ψ) && k <= length(π_nonint)
            i = i + 1
            if Ψ[i] ∈ [0, 1]
                π_nonint[r[i]] = Ψ[i]
                # replace that unit with a new unit
                Ψ[i] = π_nonint[k]
                # And also in the auxillary data matrix
                B[:, i] = A[:, k]
                # update the vector that keeps track of unit indexes
                r[i] = k

                k = k + 1
            end
        end
    end
    # now do the final iteration
    Ψ = update_psi(Ψ, B)
    π_nonint[r] = Ψ
    return (π_nonint)
end

function update_psi(Ψ, B)
    # get vector u in the kernel of B
    u = nullspace(B)[:, 1]

    # want max  λ₁, λ₂ such that 0 <= Ψ + λ₁*u <= 1 and 0 <= Ψ - λ₂*u <= 1
    # solve the inequalities for λ and you get max values for u > 0 and u < 0
    # for λ₁ : for u > 0, λ₁ = (1-Ψ)/u; for u < 0, λ₁ = -Ψ/u
    # for λ₂ : for u > 0, λ₂ = Ψ/u; for u < 0, λ₂ = (Ψ - 1)/u

    λ₁_max(; u, πₖ) = @. ifelse(u > 0, (1 - πₖ) / u, -πₖ / u)
    λ₂_max(; u, πₖ) = @. ifelse(u > 0, πₖ / u, (πₖ - 1) / u)

    λ₁_vec = filter(x -> isfinite(x), λ₁_max(; u = u, πₖ = Ψ))
    λ₂_vec = filter(x -> isfinite(x), λ₂_max(; u = u, πₖ = Ψ))

    deleteat!(λ₁_vec, λ₁_vec .<= 0)
    deleteat!(λ₂_vec, λ₂_vec .<= 0)

    λ₁ = minimum(λ₁_vec)
    λ₂ = minimum(λ₂_vec)

    # calculate the inequality expression for both lambdas
    λ₁_ineq = @. Ψ + (λ₁ * u)
    λ₂_ineq = @. Ψ - (λ₂ * u)

    # checking for floating point issues
    tol = 1e-13 # TODO: make this a sampler param
    λ₁_ineq[abs.(λ₁_ineq) .< tol] .= 0
    λ₂_ineq[abs.(λ₂_ineq) .< tol] .= 0

    # the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
    q₁ = λ₂ / (λ₁ + λ₂)

    new_πₖ = rand() < q₁ ? λ₁_ineq : λ₂_ineq
    return (new_πₖ)
end

function cubeland(pikstar, πₖ, x)
    ### Landing Phase ###
    # Goal: Find sample s such that E(s|π*) = π*, where π* is output from flight phase
    # q non-integer elements of π should be <= p auxillary variables

    pikstar_return = copy(pikstar)
    # get all non-integer probabilities
    non_int_ind = findall(u -> u .∉ Ref(Set([0, 1])), pikstar_return)
    non_int_piks = pikstar_return[non_int_ind]
    N_land = length(non_int_piks)

    # get auxillary variables for those units
    x_land = x[:, non_int_ind]

    # Get all possible samples combinations for the non-integer units
    # first, get the sample size from the total inclusion probability
    total_prob = sum(non_int_piks)
    n_land = round(Int, total_prob)

    # rescale so that the inclusion probabilities sum to the sample size
    non_int_piks = n_land * (non_int_piks / sum(non_int_piks))

    # then get matrix of potential sample design
    # get vector with appropriate allocation of 0's and 1's
    base_vec = vcat(repeat([1.0]; outer = n_land), repeat([0.0]; outer = (N_land - n_land)))

    samps = reduce(vcat, transpose.(unique_permutations(base_vec)))

    #Let's calculate the cost for each potential sampling design
    # This is C_2(s) from the appendix of Deville and Tillé 2004
    # C(s) = (s - π*)'A'(AA')^-1 A(s - π*)

    # get matrix of (s - π*), samps has a sample for each row
    sub_mat = samps .- reshape(non_int_piks, :, N_land)

    # let's get A for the non-integer units
    A_land = x_land ./ reshape(πₖ[non_int_ind], :, N_land)

    sample_pt = A_land * transpose(sub_mat)
    ## FIXME: need to deal with the case that there are fixed zeros in πₖ
    #A = x ./ reshape(πₖ, :, N)
    zero_pik_ind = findall(isequal(0), πₖ)
    A =
        x[:, setdiff(1:end, zero_pik_ind)] ./
        reshape(πₖ[setdiff(1:end, zero_pik_ind)], :, length(πₖ) - length(zero_pik_ind))

    cost = zeros(size(samps, 1))
    for i in 1:size(samps)[1]
        cost[i] = transpose(sample_pt[:, i]) * inv(A * transpose(A)) * sample_pt[:, i]
    end

    # get matrix of samples and costs
    id = 1:size(samps)[1]
    lp_mat = [id cost samps]

    ## linear programing ##
    model = Model(HiGHS.Optimizer)

    @variable(model, ps[1:size(samps, 1)] >= 0)

    # multiple cost (lp_mat[2]) by ps[id], where id is lp_mat[1]
    @objective(
        model,
        Min,
        sum(sample[2] * ps[trunc(Int, sample[1])] for sample in eachrow(lp_mat))
    )

    @constraint(model, sum(ps[id]) == 1)

    for i in 1:size(samps, 2)
        @constraint(model, sum(ps .* (samps .> 0)[:, i]) == non_int_piks[i])
    end

    optimize!(model)

    if has_values(model)
        samp_prob = value.(ps)

        # pick a sample based on their probabilities
        samp_ind = StatsBase.sample(1:length(samp_prob), Weights(samp_prob))

        # fill in non-integer points with the sample option picked by lp
        pikstar_return[non_int_ind] = samps[samp_ind, :]
    else
        @warn "The linear program did not find a feasible solution."
        pikstar_return[non_int_ind] = samps[sample(1:size(samps, 1)), :]
    end

    return pikstar_return
end

# all credit to stackoverflow https://stackoverflow.com/questions/65051953/julia-generate-all-non-repeating-permutations-in-set-with-duplicates
function unique_permutations(x::T, prefix = T()) where {T}
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i - 1]
                continue
            end
            append!(
                t,
                unique_permutations([x[begin:(i - 1)]; x[(i + 1):end]], [prefix; x[i]]),
            )
        end
        return t
    end
end

#=
"""
    CubeSampling

A `BONRefiner` that uses Cube Sampling (Tillé 2011)

**numsites**, an Integer (def. 50), specifying the number of points to use.

**fast**, a Boolean (def. true) indicating whether to use the fast flight algorithm. 

**x**, a Matrix of auxillary variables for the candidate points, with one row for each variable and one column for each candidate point.

**πₖ**, a Float Vector indicating the probabilities of inclusion for each candidate point; should sum to numsites value.
"""

Base.@kwdef struct CubeSampling{I <: Integer, M <: Matrix, V <: Vector} <: BONSampler
    numsites::I = 50
    fast::Bool = true
    x::M = rand(0:4, 3, 50)
    πₖ::V = zeros(size(x, 2))
    function CubeSampling(numsites, fast, x, πₖ)
        cs = new{typeof(numsites), typeof(x), typeof(πₖ)}(numsites, fast, x, πₖ)
        _check_arguments(cs)
        return cs
    end
end

numsites(cubesampling::CubeSampling) = cubesampling.numsites
maxsites(cubesampling::CubeSampling) = size(cubesampling.x, 2)

function check_arguments(cubesampling::CubeSampling)
    check(TooFewSites, cubesampling)
    check(TooManySites, cubesampling)

    if numsites > length(cubesampling.πₖ)
        throw(
            ArgumentError(
                "You cannot select more points than the number of candidate points.",
            ),
        )
    end
    if length(cubesampling.πₖ) != size(cubesampling.x, 2)
        throw(
            DimensionMismatch(
                "The number of inclusion probabilites does not match the dimensions of the auxillary variable matrix.",
            ),
        )
    end
    return
end

function check_conditions(coords, pool, sampler)
    πₖ = sampler.πₖ
    if sum(sampler.πₖ) == 0
        @info "Probabilities of inclusion were not provided, so we assume equal probability design."
        πₖ = [sampler.numsites / length(pool) for _ in eachindex(pool)]
    end
    if round(Int, sum(πₖ)) != sampler.numsites
        @warn "The inclusion probabilities sum to $(round(Int, sum(πₖ))), which will be your sample size instead of numsites."
    end
    if length(pool) != length(πₖ)
        throw(
            DimensionMismatch(
                "The πₖ vector does not match the number of candidate points.",
            ),
        )
    end
    if length(πₖ) != size(sampler.x, 2)
        throw(
            DimensionMismatch(
                "There is a mismatch in the number of inclusion probabilities and the points in the auxillary matrix.",
            ),
        )
    end
    return πₖ
end

function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::CubeSampling,
    uncertainty::Matrix{T},
) where {T <: AbstractFloat}
    πₖ = check_conditions(coords, pool, sampler)

    # sort points by distance in auxillary variable space
    dist = mahalanobis(πₖ, sampler.x)
    perm = sortperm(dist; rev = true)
    pool, πₖ = pool[perm], πₖ[perm]

    x = sampler.x[:, perm]

    # To keep the sample size enforced, add πₖ as an aux variable
    x = vcat(transpose(πₖ), x)

    # pick flight phase algorithm
    π_optimal_flight = sampler.fast ? cubefastflight(πₖ, x) : cubeflight(πₖ, x)
    # check if there are non-integer probabilities
    non_int_ind = findall(u -> u .∉ Ref(Set([0, 1])), π_optimal_flight)
    # if so, perform landing phase to resolve them
    π_optimal = isempty(non_int_ind) ? π_optimal_flight : cubeland(π_optimal_flight, πₖ, x)

    selected = pool[findall(isequal(1), π_optimal)]

    for i in eachindex(selected)
        coords[i] = pool[i]
    end

    return (coords, uncertainty)
end

function cubeflight(πₖ, x)
    N = length(πₖ)
    rowdim = size(x)[1]

    j = 0
    set_nullspace = zeros(1, 2)
    π_optimal = πₖ
    # check if there is a possible u to satisfy the conditions
    while size(set_nullspace)[2] != 0
        j = j + 1

        ## STEP 1 ##

        # find a vector u that is in the kernel of the matrix A
        # A is the matrix of auxillary variables didvided by the inclusion probability
        # for the population unit
        A = similar(x, Float64)
        for i in 1:N
            if π_optimal[i] .∈ Ref(Set([0, 1]))
                A[:, i] = zeros(rowdim)
            else
                A[:, i] = x[:, i] ./ π_optimal[i]
            end
        end

        # get the nullspace of A
        kernel = nullspace(A)

        # u is in the kernel of A, but also u_k = 0 when π_k is {0,1}
        # let's make sure the rows that need it satisfy that condition

        # get index where π_optimal is 0 or 1
        π_integer = findall(u -> u .∈ Ref(Set([0, 1])), π_optimal)

        # if none of the π_optimal's are fixed yet (as 0 or 1) u can be a vector from the nullspace
        if length(π_integer) == 0
            u = kernel[:, rand(1:size(kernel)[2])]

            # if only one is fixed, can also pick a u vector but it shouldn't be the trivial indicator vector
        elseif length(π_integer) == 1
            sums = sum(eachrow(kernel))
            # find indicator vector
            ind = findall(isequal(1), sums)

            # get vector of potential column indices, remove unit column, and get random u
            ind = deleteat!(collect(1:size(kernel)[2]), ind)
            u = kernel[:, rand(ind)]

            # otherwise, need to make sure u_k = 0 condition is satisfied for fixed pikstar's
        else
            # get rows of A's nullspace corresponding to those pikstar's
            set_A = kernel[π_integer, :]
            # get the nullspace of that matrix

            set_nullspace = nullspace(set_A)

            if size(set_nullspace)[2] == 0
                break
            end

            # randomly pick a vector from the second nullspace
            v = set_nullspace[:, rand(1:size(set_nullspace)[2])]

            # multiply it by original kernel to get a vector u that satisfies the zeroes requirement
            ### FIX ME: need to deal with rounding issues, the zeros are not real zeroes!
            u = kernel * v
            # this is a hacky way to make sure zeros are real zeroes 
            u[π_integer] .= 0
        end

        ## STEP 2 ##

        # Find the maximum λ₁, λ₂ such that:
        # 1) -πₖ⁺ <= λ₁*u <= 1 - πₖ⁺ 
        # 2) -πₖ⁺ <= -λ₂*u <= 1 - πₖ⁺
        # 
        # Solve the inequalities for λ and you get max values for cases where u > 0 and u < 0
        # for u > 0: 
        #   λ₁ = (1-πₖ⁺)/u; 
        #   λ₂ = πₖ⁺/u; 
        # for u < 0: 
        #    λ₁ = -πₖ⁺/u
        #    λ₂ = (πₖ⁺ - 1)/u

        λ₁_max(; u, πₖ) = @. ifelse(u > 0, (1 - πₖ) / u, -πₖ / u)
        λ₂_max(; u, πₖ) = @. ifelse(u > 0, πₖ / u, (πₖ - 1) / u)

        λ₁ = minimum(filter(isfinite, λ₁_max(; u = u, πₖ = πₖ)))
        λ₂ = minimum(filter(isfinite, λ₂_max(; u = u, πₖ = πₖ)))

        ## STEP 3 ##

        # calculate the inequality expression for both lambdas
        λ1_ineq = @. πₖ + (λ₁ * u)
        λ2_ineq = @. πₖ - (λ₂ * u)

        ineq_mat = reduce(hcat, [λ1_ineq, λ2_ineq])
        # the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
        q₁ = λ₂ / (λ₁ + λ₂)
        q₂ = 1 - q₁

        π_optimal = map(r -> sample(r, Weights([q₁, q₂])), eachrow(ineq_mat))
    end
    return (π_optimal)
end

function cubefastflight(πₖ, x)

    # number of auxillary variables
    num_aux = size(x)[1]

    # get all non-integer probabilities
    non_int_ind = findall(u -> u .∉ Ref(Set([0, 1])), πₖ)

    π_nonint = πₖ[non_int_ind]

    Ψ = π_nonint[1:(num_aux + 1)]
    r = collect(1:(num_aux + 1))

    A = x[:, non_int_ind] ./ reshape(π_nonint, :, length(non_int_ind))
    B = A[:, 1:(num_aux + 1)]

    k = num_aux + 2

    while k <= length(π_nonint)
        Ψ = update_psi(Ψ, B)

        if length(findall(z -> z .< 0, Ψ)) > 0
            throw(error("Negative inclusion probability"))
        end

        # update for the probabilities that are now integers
        i = 0
        while i < length(Ψ) && k <= length(π_nonint)
            i = i + 1
            if Ψ[i] ∈ [0, 1]
                π_nonint[r[i]] = Ψ[i]
                # replace that unit with a new unit
                Ψ[i] = π_nonint[k]
                # And also in the auxillary data matrix
                B[:, i] = A[:, k]
                # update the vector that keeps track of unit indexes
                r[i] = k

                k = k + 1
            end
        end
    end
    # now do the final iteration
    Ψ = update_psi(Ψ, B)
    π_nonint[r] = Ψ
    return (π_nonint)
end

function update_psi(Ψ, B)
    # get vector u in the kernel of B
    u = nullspace(B)[:, 1]

    # want max  λ₁, λ₂ such that 0 <= Ψ + λ₁*u <= 1 and 0 <= Ψ - λ₂*u <= 1
    # solve the inequalities for λ and you get max values for u > 0 and u < 0
    # for λ₁ : for u > 0, λ₁ = (1-Ψ)/u; for u < 0, λ₁ = -Ψ/u
    # for λ₂ : for u > 0, λ₂ = Ψ/u; for u < 0, λ₂ = (Ψ - 1)/u

    λ₁_max(; u, πₖ) = @. ifelse(u > 0, (1 - πₖ) / u, -πₖ / u)
    λ₂_max(; u, πₖ) = @. ifelse(u > 0, πₖ / u, (πₖ - 1) / u)

    λ₁_vec = filter(x -> isfinite(x), λ₁_max(; u = u, πₖ = Ψ))
    λ₂_vec = filter(x -> isfinite(x), λ₂_max(; u = u, πₖ = Ψ))

    deleteat!(λ₁_vec, λ₁_vec .<= 0)
    deleteat!(λ₂_vec, λ₂_vec .<= 0)

    λ₁ = minimum(λ₁_vec)
    λ₂ = minimum(λ₂_vec)

    # calculate the inequality expression for both lambdas
    λ₁_ineq = @. Ψ + (λ₁ * u)
    λ₂_ineq = @. Ψ - (λ₂ * u)

    # checking for floating point issues
    tol = 1e-13
    λ₁_ineq[abs.(λ₁_ineq) .< tol] .= 0
    λ₂_ineq[abs.(λ₂_ineq) .< tol] .= 0

    # the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
    q₁ = λ₂ / (λ₁ + λ₂)

    new_πₖ = rand() < q₁ ? λ₁_ineq : λ₂_ineq
    return (new_πₖ)
end

function cubeland(pikstar, πₖ, x)
    ### Landing Phase ###
    # Goal: Find sample s such that E(s|π*) = π*, where π* is output from flight phase
    # q non-integer elements of π should be <= p auxillary variables

    pikstar_return = copy(pikstar)
    # get all non-integer probabilities
    non_int_ind = findall(u -> u .∉ Ref(Set([0, 1])), pikstar_return)
    non_int_piks = pikstar_return[non_int_ind]
    N_land = length(non_int_piks)

    # get auxillary variables for those units
    x_land = x[:, non_int_ind]

    # Get all possible samples combinations for the non-integer units
    # first, get the sample size from the total inclusion probability
    total_prob = sum(non_int_piks)
    n_land = round(Int, total_prob)

    # rescale so that the inclusion probabilities sum to the sample size
    non_int_piks = n_land * (non_int_piks / sum(non_int_piks))

    # then get matrix of potential sample design
    # get vector with appropriate allocation of 0's and 1's
    base_vec = vcat(repeat([1.0]; outer = n_land), repeat([0.0]; outer = (N_land - n_land)))

    samps = reduce(vcat, transpose.(unique_permutations(base_vec)))

    #Let's calculate the cost for each potential sampling design
    # This is C_2(s) from the appendix of Deville and Tillé 2004
    # C(s) = (s - π*)'A'(AA')^-1 A(s - π*)

    # get matrix of (s - π*), samps has a sample for each row
    sub_mat = samps .- reshape(non_int_piks, :, N_land)

    # let's get A for the non-integer units
    A_land = x_land ./ reshape(πₖ[non_int_ind], :, N_land)

    sample_pt = A_land * transpose(sub_mat)
    ## FIXME: need to deal with the case that there are fixed zeros in πₖ
    #A = x ./ reshape(πₖ, :, N)
    zero_pik_ind = findall(isequal(0), πₖ)
    A =
        x[:, setdiff(1:end, zero_pik_ind)] ./
        reshape(πₖ[setdiff(1:end, zero_pik_ind)], :, length(πₖ) - length(zero_pik_ind))

    cost = zeros(size(samps, 1))
    for i in 1:size(samps)[1]
        cost[i] = transpose(sample_pt[:, i]) * inv(A * transpose(A)) * sample_pt[:, i]
    end

    # get matrix of samples and costs
    id = 1:size(samps)[1]
    lp_mat = [id cost samps]

    ## linear programing ##
    model = Model(HiGHS.Optimizer)

    @variable(model, ps[1:size(samps, 1)] >= 0)

    # multiple cost (lp_mat[2]) by ps[id], where id is lp_mat[1]
    @objective(
        model,
        Min,
        sum(sample[2] * ps[trunc(Int, sample[1])] for sample in eachrow(lp_mat))
    )

    @constraint(model, sum(ps[id]) == 1)

    for i in 1:size(samps, 2)
        @constraint(model, sum(ps .* (samps .> 0)[:, i]) == non_int_piks[i])
    end

    optimize!(model)

    if has_values(model)
        samp_prob = value.(ps)

        # pick a sample based on their probabilities
        samp_ind = sample(1:length(samp_prob), Weights(samp_prob))

        # fill in non-integer points with the sample option picked by lp
        pikstar_return[non_int_ind] = samps[samp_ind, :]
    else
        @warn "The linear program did not find a feasible solution."
        pikstar_return[non_int_ind] = samps[sample(1:size(samps, 1)), :]
    end

    return pikstar_return
end

# all credit to stackoverflow https://stackoverflow.com/questions/65051953/julia-generate-all-non-repeating-permutations-in-set-with-duplicates
function unique_permutations(x::T, prefix = T()) where {T}
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i - 1]
                continue
            end
            append!(
                t,
                unique_permutations([x[begin:(i - 1)]; x[(i + 1):end]], [prefix; x[i]]),
            )
        end
        return t
    end
end

function mahalanobis(πₖ, x)
    # drop variables that are the same for all points
    num_uniq = map(y -> length(unique(y)), eachrow(x))
    nonuniq_ind = findall(z -> z == 1, num_uniq)
    x = length(nonuniq_ind) > 0 ? x[1:end .!= nonuniq_ind, :] : x

    N = length(πₖ)
    p = size(x, 1)

    x̂ = x ./ reshape(πₖ, :, N)
    mean_x̂ = (1 / N) .* sum(x̂; dims = 2)

    k_vecs = x̂ .- mean_x̂
    outer_prods = Array{Float64}(undef, p, p, N)
    for i in 1:N
        outer_prods[:, :, i] = k_vecs[1:p, i] * transpose(k_vecs[1:p, i])
    end

    sigma = (1 / (N - 1)) * dropdims(sum(outer_prods; dims = 3); dims = 3)
    inv_sigma = inv(sigma)

    d = Vector{Float64}(undef, N)
    for i in 1:N
        d[i] = (transpose(x̂[1:p, i] - mean_x̂) * inv_sigma * (x̂[1:p, i] - mean_x̂))[1]
    end

    return d
end 

=#