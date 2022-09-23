using LinearAlgebra
using StatsBase

# Need a vector of inclusion probabilities for a set of points


# Let's implement a little concrete example using the notation from Tillé 2011
# I use the same structure as the paper, with two phases Flight and Landing 

# n samples in sample s
n = 5
# from population U of size N
N = 10
# with probability of inclusion
pik = rand(Float64, 10)

# p auxillary variables
p = 3

# generate random auxillary variables, p rows are x_k variables for N sample units
x = rand(0:4, p, N)

### Flight Phase ###
 i = 0
 set_nullspace = zeros(1,2)
# check if there is a possible u to satisfy the conditions
while size(set_nullspace)[2] != 0
    i = i+1
    #yield(i)
    ## STEP 1 ##

    # find a vector u that is in the kernel of the matrix A
    # A is the matrix of auxillary variables didvided by the inclusion probability
    # for the population unit
    A = similar(x, Float64)
    for i = 1:N
        if pik[i] .∈ Ref(Set([0,1]))
            A[:,i] = zeros(p) # p is the row dimension
        else 
            A[:,i] = x[:,i] ./ pik[i]
        end
    end

    # get the nullspace of A
    kernal = nullspace(A)

    # u is in the kernal of A, but also u_k = 0 when π_k is {0,1}
    # let's make sure the rows that need it satisfy that condition

    # get index where pik is 0 or 1
    set_piks = findall(x -> x .∈ Ref(Set([0,1])), pik)

    # if none of the pik's are fixed yet (as 0 or 1) u can be a vector from the nullspace
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

    # otherwise, need to make sure u_k = 0 condition is satisfied for fixed pik's
    else
        # get rows of A's nullspace corresponding to those pik's
        set_A = kernal[set_piks, :]
        # get the nullspace of that matrix
        print(i)
        set_nullspace = nullspace(set_A)
        print(i)

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

    # want max  λ_1, λ_2 such that -pik <= λ_1 * u <= 1 - pik and -pik <= -λ_2 * u <= 1 - pik
    # solve the inequalities for λ and you get max values for u > 0 and u < 0
    # for λ_1 : for u > 0, λ_1 = (1-pik)/u; for u < 0, λ_1 = -pik/u
    # for λ_2 : for u > 0, λ_2 = pik/u; for u < 0, λ_2 = (pik - 1)/u

    λ1_max(; u, pik) = @. ifelse(u > 0, (1 - pik) / u, - pik / u)
    λ2_max(; u, pik) = @. ifelse(u > 0, pik / u, (pik - 1) / u)

    #vars_df = DataFrame(pik = pik, u = u, λ1 = λ1_max(u = u, pik = pik), λ2 = λ2_max(u = u, pik = pik)) 
    λ1 = minimum(filter(x -> isfinite(x), λ1_max(u = u, pik = pik)))
    λ2 = minimum(filter(x -> isfinite(x), λ2_max(u = u, pik = pik)))

    ## STEP 3 ##

    # calculate the inequality expression for both lambdas
    λ1_ineq = @. pik + ( λ1 * u)
    λ2_ineq = @. pik - ( λ2 * u)

    ineq_mat = reduce(hcat, [λ1_ineq, λ2_ineq])
    # the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
    q1 = λ2/(λ1 + λ2)
    q2 = 1 - q1

    pik = map(r -> sample(r, Weights([q1,q2])), eachrow(ineq_mat))
    #yield(i)

end 
# double check that it satisfies the inequality 
#= if minimum(λ1_ineq) < 0
    print("λ1 lower bound problem") 
end
if maximum(λ1_ineq)  > 1
    print("λ1 upper bound problem") 
end
if minimum(λ2_ineq)  < 0
    print("λ2 lower bound problem") 
end
if maximum(λ2_ineq) > 1
    print("λ2 upper bound problem") 
end =#

