using LinearAlgebra
using DataFrames
using Pipe

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

## STEP 1 ##

# find a vector u that is in the kernel of the matrix A
# A is the matrix of auxillary variables didvided by the inclusion probability
# for the population unit
A = similar(x, Float64)
for i = 1:N
    A[:,i] = x[:,i] ./ pik[i]
end

## FIXME: this does not deal with the case of a zero or one inclusion probability
##        which should result in u_k(t) = 0, and currenlty will divide by zero to get A

# get the nullspace of A
kernal = nullspace(A)

# randomly pick a vector u in that kernel that is non-zero
# get vector of column sums
sums = sum(eachrow(kernal))
# find zeros
zero_ind = findall(x -> x==0, sums)

# get vector of potential column indexes, remove zero columns, and get random u
ind = deleteat!(collect(1:N), zero_ind)
u = kernal[:, rand(ind)]

## STEP 2 ##

# want max  λ_1, λ_2 such that -pik <= λ_1 * u <= 1 - pik and -pik <= -λ_2 * u <= 1 - pik
# solve the inequalities for λ and you get max values for u > 0 and u < 0
# for λ_1 : for u > 0, λ_1 = (1-pik)/u; for u < 0, λ_1 = -pik/u
# for λ_2 : for u > 0, λ_2 = pik/u; for u < 0, λ_2 = (pik - 1)/u

λ1_max(; u, pik) = @. ifelse(u > 0, (1 - pik) / u, - pik / u)
λ2_max(; u, pik) = @. ifelse(u > 0, pik / u, (pik - 1) / u)

#vars_df = DataFrame(pik = pik, u = u, λ1 = λ1_max(u = u, pik = pik), λ2 = λ2_max(u = u, pik = pik)) 
λ1 = minimum(λ1_max(u = u, pik = pik))
λ2 = minimum(λ2_max(u = u, pik = pik))

## STEP 3 ##

# calculate the inequality expression for both lambdas
λ1_ineq = @. pik + ( λ1 * u)
λ2_ineq = @. pik - ( λ2 * u)

ineq_mat = reduce(hcat, [λ1_ineq, λ2_ineq])
# the new inclusion probability π is one of the lambda expressions with a given probability q1, q2
q1 = λ2/(λ1 + λ2)
q2 = 1 - q1

new_pik = sample.(eachrow(ineq_mat), [Weights([q1,q2])])

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

