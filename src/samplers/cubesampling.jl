"""
    CubeSampling

`CubeSampling` implements the cube method of [Deville2004EffSam](@cite) for
balanced sampling with respect to a set of auxiliary variables (features).

The algorithm proceeds in two phases:
- Flight phase: probabilities are iteratively moved towards 0/1 while preserving
  linear constraints on the mean of auxiliary variables in expectation (including fixed sample size).
- Landing phase: if fractional probabilities remain, an optimization step chooses
  a 0/1 sample that best matches the target constraints.

If `inclusion` is not provided, uniform inclusion probabilities are derived from
`sampler.num_nodes` and the domain pool size. Returned nodes correspond to units
with final probability 1.
"""
@kwdef struct CubeSampling <: BONSampler
    num_nodes = _DEFAULT_NUM_NODES
end

"""
    _sample(sampler::CubeSampling, domain; inclusion=nothing)

Draw a sample using the cube method, balancing means of auxiliary variables
(`getfeatures(domain)`) while achieving the desired sample size in expectation.

Arguments:
- `sampler.num_nodes`: desired number of selected sites
- `domain`: sampling domain; must support `getpool(domain)` and `getfeatures(domain)`
- `inclusion`: optional vector/array of inclusion probabilities indexed by pool items

Returns a `BiodiversityObservationNetwork`.
"""
function _sample(
    sampler::CubeSampling,
    domain;
    inclusion = nothing
)
    inclusion = isnothing(inclusion) ? get_uniform_inclusion(sampler, domain) : inclusion

    inclusion = [inclusion[p] for p in getpool(domain)]

    features = getfeatures(domain)

    # Sort domain units by Mahalanobis distance to mean feature vector
    sorted_idx = _sort_features_by_mahalanobis(features, inclusion)

    sorted_inclusion = inclusion[sorted_idx]
    sorted_features = features[:, sorted_idx]

    # Flight phase of cube method
    π_flight = _cube_flight_phase(sorted_inclusion, sorted_features)
    
    # Landing phase if fractional probabilities remain
    π_optimal = any(x -> x ∉ [0,1], π_flight) ? _cube_landing_phase(π_flight, sorted_inclusion, sorted_features) : π_flight

    # Select included nodes from sorted pool
    sorted_pool = getpool(domain)[sorted_idx]
    
    selected_idx = Bool.(π_optimal)
    nodes = sorted_pool[selected_idx]
    selected_features = features[:,selected_idx]

    return BiodiversityObservationNetwork(nodes, selected_features)
end 


"""
    _sort_features_by_mahalanobis(features, inclusion)

Order units to stabilize the flight phase by spreading early decisions across
feature space. Units are sorted by Mahalanobis distance from the inclusion-
weighted mean feature vector.
"""
function _sort_features_by_mahalanobis(features, inclusion)
    # Mahalanobis distance matrix between feature vectors
    mh = Distances.Mahalanobis(features * features')
    
    # Normalize features by inclusion probabilities
    x = features ./ inclusion'
    x_bar = SB.mean(eachcol(x))

    # Sort features by Mahalanobis distance to mean (descending)
    perm = sortperm([mh(x_bar, x) for x in eachcol(features)]) |> reverse
    return perm
end


"""
    _cube_flight_phase(πₖ, x)

Run the flight phase of the cube method.

`πₖ` are current inclusion probabilities; `x` (auxiliary matrix) is augmented
with a first row of `πₖ'` so sample size fixed. The method repeatedly finds a direction in the null space of the constraint matrix and pushes a small subset of probabilities to 0 or 1 while preserving
balances in expectation.
"""
function _cube_flight_phase(πₖ, x)
    x = vcat(transpose(πₖ), x) # add inclusion probabilities to fix sample size

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
        Ψ =_update_psi(Ψ, B)

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
    Ψ = _update_psi(Ψ, B)
    π_nonint[r] = Ψ
    return π_nonint
end

"""
    _update_psi(Ψ, B)

Given current working vector `Ψ` and constraint block `B`, compute a feasible
update along a null-space direction `u` by maximizing step sizes `λ₁, λ₂` that
keep probabilities within [0,1].
"""
function _update_psi(Ψ, B)
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

"""
    _cube_landing_phase(pikstar, πₖ, x)

Run the landing phase if some probabilities are still fractional after flight.
Formulate a small linear program over the fractional units to select a 0/1
sample whose auxiliary totals match `pikstar` in expectation with minimal cost
`C(s)` (per Deville & Tillé, 2004).
"""
function _cube_landing_phase(pikstar, πₖ, x)
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

    # Let's calculate the cost for each potential sampling design
    # This is C_2(s) from the appendix of Deville and Tillé 2004
    # C(s) = (s - π*)'A'(AA')^-1 A(s - π*)

    # get matrix of (s - π*), samps has a sample for each row
    sub_mat = samps .- reshape(non_int_piks, :, N_land)

    # let's get A for the non-integer units
    A_land = x_land ./ reshape(πₖ[non_int_ind], :, N_land)

    sample_pt = A_land * transpose(sub_mat)

    ## Deal with the case that there are fixed zeros in πₖ
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

"""
    unique_permutations(x::T, prefix = T()) where {T}

Generate all unique permutations for a multiset `x` without repetition of
duplicates. Based on StackOverflow
(`https://stackoverflow.com/questions/65051953/julia-generate-all-non-repeating-permutations-in-set-with-duplicates`).
"""
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
