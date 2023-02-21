"""
    CubeSampling

A `BONRefiner` that uses Cube Sampling (Tillé 2011)

...

**numpoints**, an Integer (def. 50), specifying the number of points to use.

**fast**, a Boolean (def. true) indicating whether to use the fast flight algorithm. 

**x**, a Matrix of auxillary variables for the candidate points, with one row for each variable and one column for each candidate point.

**pik**, a Float Vector indicating the probabilities of inclusion for each candidate point; should sum to numpoints value.
"""

Base.@kwdef mutable struct CubeSampling{I <: Integer, M <: Matrix, V <: Vector} <: BONRefiner
    numpoints::I = 50
    fast::Bool = true
    x::M = rand(0:4, 3, 50)
    pik::V = zeros(size(x,2))

    function CubeSampling(numpoints, fast, x, pik)

        if numpoints < one(numpoints)
            throw(ArgumentError("You cannot have a CubeSampling with fewer than one point.",),)
        end
        if numpoints > length(pik)
            throw(ArgumentError("You cannot select more points than the number of candidate points.",),)
        end
        if length(pik) != size(x, 2)
            throw(DimensionMismatch("The number of inclusion probabilites does not match the dimensions of the auxillary variable matrix.",),)
        end
        return new{typeof(numpoints), typeof(x), typeof(pik)}(numpoints, fast, x, pik)
    end
end

function _generate!(
    coords::Vector{CartesianIndex}, 
    pool::Vector{CartesianIndex},
    sampler::CubeSampling,
    uncertainty::Matrix{T}
    ) where {T <: AbstractFloat}
    
    #check if they gave us pik or not
    if sum(sampler.pik) == 0
        @info "Probabilities of inclusion were not provided, so we assume equal probability design."
        pik = fill(sampler.numpoints/length(pool), length(pool))
    else pik = sampler.pik
    end

    if sum(pik) != sampler.numpoints
        throw(Warning("The inclusion probabilities sum to $sum(pik), which will be your sample size instead of numpoints."))
    # check that dimensions match
    if length(pool) != length(pik)
        throw(DimensionMismatch("The pik vector does not match the number of candidate points.",),)
    end

    if length(pik) != size(sampler.x, 2)
        throw(DimensionMismatch("There is a mismatch in the number of inclusion probabilities and the points in the auxillary matrix.",),)
    end

    # sort points by distance in auxillary variable space
    dist = mahalanobis(pik, sampler.x)
    perm = sortperm(dist, rev=true)

    pool = pool[perm]
    pik = pik[perm]
    x = sampler.x[:,perm]
 
    # if we want the sample size enforced, add pik as an aux variable
    x = vcat(transpose(pik), x)

    # pick flight phase algorithm
    pikstar = sampler.fast ? cubefastflight(pik, x) : cubeflight(pik, x)
    # check if there are non-integer probabilities
    non_int_ind = findall(x -> x .∉ Ref(Set([0,1])), pikstar)
    # if so, perform landing phase to resolve them
    pikstar = isempty(non_int_ind) ? pikstar : cubeland(pikstar, pik, x)

    selected = pool[findall(x -> x == 1, pikstar)]

    for i = 1:length(selected)
        coords[i] = pool[i]
    end

    return (coords, uncertainty)
end

function cubeflight(pik, x)
    
    N = length(pik)
    n = sum(pik)
    p = size(x)[1]

    ### Flight Phase ###
    j = 0
    set_nullspace = zeros(1,2)
    pikstar = pik
    # check if there is a possible u to satisfy the conditions
    while size(set_nullspace)[2] != 0
        j = j+1

        ## STEP 1 ##

        # find a vector u that is in the kernel of the matrix A
        # A is the matrix of auxillary variables didvided by the inclusion probability
        # for the population unit
        A = similar(x, Float64)
        for i = 1:N
            if pikstar[i] .∈ Ref(Set([0,1]))
                A[:,i] = zeros(p) # p is the row dimension
            else 
                A[:,i] = x[:,i] ./ pikstar[i]
            end
        end

        # get the nullspace of A
        kernal = nullspace(A)

        # u is in the kernal of A, but also u_k = 0 when π_k is {0,1}
        # let's make sure the rows that need it satisfy that condition

        # get index where pikstar is 0 or 1
        set_piks = findall(x -> x .∈ Ref(Set([0,1])), pikstar)

        # if none of the pikstar's are fixed yet (as 0 or 1) u can be a vector from the nullspace
        if length(set_piks) == 0
            u = kernal[:, rand(1:size(kernal)[2])]
        
        # if only one is fixed, can also pick a u vector but it shouldn't be the trivial indicator vector
        elseif length(set_piks) == 1
            sums = sum(eachrow(kernal))
            # find indicator vector
            ind = findall(x -> x==1, sums)

            # get vector of potential column indices, remove unit column, and get random u
            ind = deleteat!(collect(1:size(kernal)[2]), ind)
            u = kernal[:, rand(ind)]

        # otherwise, need to make sure u_k = 0 condition is satisfied for fixed pikstar's
        else
            # get rows of A's nullspace corresponding to those pikstar's
            set_A = kernal[set_piks, :]
            # get the nullspace of that matrix
            
            set_nullspace = nullspace(set_A)

            if size(set_nullspace)[2] == 0
                break
            end

            # randomly pick a vector from the second nullspace
            v = set_nullspace[:, rand(1:size(set_nullspace)[2])]

            # multiply it by original kernal to get a vector u that satisfies the zeroes requirement
            ### FIX ME: need to deal with rounding issues, the zeros are not real zeroes!
            u = kernal * v
            # this is a hacky way to make sure zeros are real zeroes 
            u[set_piks] .= 0
        end

        ## STEP 2 ##

        # want max  λ_1, λ_2 such that -pikstar <= λ_1 * u <= 1 - pikstar and -pikstar <= -λ_2 * u <= 1 - pikstar
        # solve the inequalities for λ and you get max values for u > 0 and u < 0
        # for λ_1 : for u > 0, λ_1 = (1-pikstar)/u; for u < 0, λ_1 = -pikstar/u
        # for λ_2 : for u > 0, λ_2 = pikstar/u; for u < 0, λ_2 = (pikstar - 1)/u

        λ1_max(; u, pik) = @. ifelse(u > 0, (1 - pik) / u, - pik / u)
        λ2_max(; u, pik) = @. ifelse(u > 0, pik / u, (pik - 1) / u)

        #vars_df = DataFrame(pikstar = pikstar, u = u, λ1 = λ1_max(u = u, pikstar = pikstar), λ2 = λ2_max(u = u, pikstar = pikstar)) 
        λ1 = minimum(filter(x -> isfinite(x), λ1_max(u = u, pik = pikstar)))
        λ2 = minimum(filter(x -> isfinite(x), λ2_max(u = u, pik = pikstar)))

        ## STEP 3 ##

        # calculate the inequality expression for both lambdas
        λ1_ineq = @. pikstar + ( λ1 * u)
        λ2_ineq = @. pikstar - ( λ2 * u)

        ineq_mat = reduce(hcat, [λ1_ineq, λ2_ineq])
        # the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
        q1 = λ2/(λ1 + λ2)
        q2 = 1 - q1

        pikstar = map(r -> sample(r, Weights([q1,q2])), eachrow(ineq_mat))
        #yield(i)

    end 
    return(pikstar)
end

function cubefastflight(pik, x)

    ## Initialization ##

    # number of auxillary variables
    p = size(x)[1]

    # get all non-integer probabilities
    non_int_ind = findall(x -> x .∉ Ref(Set([0, 1])), pik)
    π = pik[non_int_ind]
    Ψ = π[1:(p + 1)]
    r = collect(1:(p + 1))

    A = x[:, non_int_ind] ./ reshape(π, :, length(non_int_ind))
    B = A[:, 1:(p + 1)]

    k = p + 2

    while k <= length(π)

        Ψ = update_psi(Ψ, B) 

        if length(findall(x -> x .<0, Ψ)) > 0
            throw(error("Negative inclusion probability"))
        end

        # update for the probabilities that are now integers
        i = 0
        while i < length(Ψ) && k <= length(π)
            i = i + 1
            if Ψ[i] ∈ [0, 1]
                # update π
                π[r[i]] = Ψ[i]
                # replace that unit with a new unit
                Ψ[i] = π[k]
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
    π[r] = Ψ

    return(π)
end

function update_psi(Ψ, B)
     # get vector u in the kernel of B
     u = nullspace(B)[:, 1]

     # want max  λ_1, λ_2 such that 0 <= Ψ + λ_1*u <= 1 and 0 <= Ψ - λ_2*u <= 1
     # solve the inequalities for λ and you get max values for u > 0 and u < 0
     # for λ_1 : for u > 0, λ_1 = (1-Ψ)/u; for u < 0, λ_1 = -Ψ/u
     # for λ_2 : for u > 0, λ_2 = Ψ/u; for u < 0, λ_2 = (Ψ - 1)/u

     λ1_max(; u, pik) = @. ifelse(u > 0, (1 - pik) / u, -pik / u)
     λ2_max(; u, pik) = @. ifelse(u > 0, pik / u, (pik - 1) / u)

     λ1_vec = filter(x -> isfinite(x), λ1_max(; u = u, pik = Ψ))
     λ2_vec = filter(x -> isfinite(x), λ2_max(; u = u, pik = Ψ))

     deleteat!(λ1_vec, λ1_vec .<= 0)
     deleteat!(λ2_vec, λ2_vec .<= 0)

     λ1 = minimum(λ1_vec)
     λ2 = minimum(λ2_vec)

     # calculate the inequality expression for both lambdas
     λ1_ineq = @. Ψ + (λ1 * u)
     λ2_ineq = @. Ψ - (λ2 * u)

     # checking for floating point issues
     tol = 1e-13
     λ1_ineq[abs.(λ1_ineq) .< tol] .= 0
     λ2_ineq[abs.(λ2_ineq) .< tol] .= 0

     ineq_mat = reduce(hcat, [λ1_ineq, λ2_ineq])
     # the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
     q1 = λ2 / (λ1 + λ2)
     q2 = 1 - q1

     Ψ = map(r -> sample(r, Weights([q1, q2])), eachrow(ineq_mat))

     return(Ψ)
end    

function cubeland(pikstar, pik, x)
    ### Landing Phase ###
    # Goal: Find sample s such that E(s|π*) = π*, where π* is output from flight phase
    # q non-integer elements of π should be <= p auxillary variables

    # get all non-integer probabilities
    non_int_ind = findall(x -> x .∉ Ref(Set([0,1])), pikstar)
    non_int_piks = pikstar[non_int_ind]
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
    base_vec = vcat(repeat([1.0], outer = n_land), repeat([0.0], outer = (N_land - n_land)))

    samps = reduce(vcat, transpose.(unique_permutations(base_vec)))

    #Let's calculate the cost for each potential sampling design
    # This is C_2(s) from the appendix of Deville and Tillé 2004
    # C(s) = (s - π*)'A'(AA')^-1 A(s - π*)

    # get matrix of (s - π*), samps has a sample for each row
    sub_mat = samps .- reshape(non_int_piks, :, N_land)

    # let's get A for the non-integer units
    A_land = x_land ./ reshape(pik[non_int_ind], :, N_land)
    

    sample_pt = A_land * transpose(sub_mat)
    ## FIXME: need to deal with the case that there are fixed zeros in pik
    #A = x ./ reshape(pik, :, N)
    zero_pik_ind = findall(x -> x == 0, pik)
    A = x[:, setdiff(1:end, zero_pik_ind)] ./ reshape(pik[setdiff(1:end, zero_pik_ind)], :, length(pik) - length(zero_pik_ind))
    

    cost = zeros(size(samps,1))
    for i in 1:size(samps)[1]
        cost[i] = transpose(sample_pt[:, i]) * inv(A*transpose(A)) * sample_pt[:, i]
    end

    # get matrix of samples and costs
    id = 1:size(samps)[1]
    lp_mat = [id cost samps]

    ## linear programing ##
    model = Model(HiGHS.Optimizer)

    @variable(model, ps[1:size(samps,1)] >= 0)

    # multiple cost (lp_mat[2]) by ps[id], where id is lp_mat[1]
    @objective(model, Min, sum(sample[2] * ps[trunc(Int, sample[1])] for sample in eachrow(lp_mat)))

    @constraint(model, sum(ps[id]) == 1)

    for i in 1:size(samps,2)
        @constraint(model, sum(ps .* (samps.>0)[:,i]) == non_int_piks[i])
    end

    optimize!(model)
    #solution_summary(model)

    #has_values(model) || @warn "The linear program did not find a feasible solution."
    if has_values(model) 
        samp_prob = value.(ps)
        
        # pick a sample based on their probabilities
        samp_ind = sample(1:length(samp_prob), Weights(samp_prob))

        # fill in non-integer points with the sample option picked by lp
        pikstar[non_int_ind] = samps[samp_ind, :]
    else
        @warn "The linear program did not find a feasible solution."
        pikstar[non_int_ind] = samps[sample(1:size(samps,1)), :]

    end 

    return(pikstar)
end

# all credit to stackoverflow https://stackoverflow.com/questions/65051953/julia-generate-all-non-repeating-permutations-in-set-with-duplicates
function unique_permutations(x::T, prefix=T()) where T
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i-1]
                continue
            end
            append!(t, unique_permutations([x[begin:i-1];x[i+1:end]], [prefix; x[i]]))
        end
        return t
    end
end

function mahalanobis(pik, x)
    # drop variables that are the same for all points
    num_uniq = map(y -> length(unique(y)), eachrow(x))
    nonuniq_ind = findall(z -> z == 1, num_uniq)
    x = length(nonuniq_ind) > 0 ? x[1:end .!= nonuniq_ind, :] : x

    N = length(pik)
    p = size(x, 1)

    x_hat = x ./ reshape(pik, :, N)
    x_hat_bar = (1/N) .* sum(x_hat, dims = 2)

    k_vecs = x_hat .- x_hat_bar
    outer_prods = Array{Float64}(undef, p, p, N)
    for i in 1:N
        outer_prods[:, :, i] = k_vecs[1:p, i] * transpose(k_vecs[1:p, i])
    end
    
    sigma = (1/(N-1)) * dropdims(sum(outer_prods, dims = 3), dims = 3)
    inv_sigma = inv(sigma)
    
    d = Vector{Float64}(undef, N)
    for i in 1:N
        d[i] = (transpose(x_hat[1:p, i] - x_hat_bar) * inv_sigma * (x_hat[1:p, i] - x_hat_bar))[1]
    end

    return d
end