"""
    SpatiallyStratified
"""
@kwdef struct SpatiallyStratified{I<:Integer,L<:Layer,F<:AbstractFloat,N} <: BONSampler
    numsites::I = 50
    layer::L = _default_strata((50, 50))
    category_dict::Dict{I,N} = _default_categories()
    inclusion_dict::Dict{N,F} = _default_inclusion()
    function SpatiallyStratified(
        numsites::I, 
        layer::L, 
        category_dict::Dict{I,N},
        inclusion_dict::Dict{N,F}
    ) where {I<:Integer,L<:Layer,F<:AbstractFloat,N}
        ss = new{I,L,F,N}(
            numsites,
            layer,
            category_dict,
            inclusion_dict
        )
        check_arguments(ss)
        return ss
    end
end


function check_arguments(ss::SpatiallyStratified)
    check(TooFewSites, ss)

    length(ss.category_dict) == length(ss.inclusion_dict) || throw(
        ArgumentError(
            "Inclusion probability vector does not have the same number of strata as there are unique values in the strata matrix",
        ),
    )

    return sum(values(ss.inclusion_dict)) ≈ 1.0 ||
           throw(ArgumentError("Inclusion probabilities across all strata do not sum to 1."))
end


function _sample!(
    coords::S,
    candidates::C,
    sampler::SpatiallyStratified{I,L,F,N},
) where {S<:Sites,C<:Sites,I,L,F,N}
    n = sampler.numsites
    strata = sampler.layer
    categories = sampler.category_dict
    
    category_ids = sort(collect(keys(categories)))
    candidate_ids = [strata.layer[x] for x in coordinates(candidates)]

    cat_idx = Dict()
    inclusion_probs = F[]
    for k in category_ids
        cat_idx[k] = findall(isequal(k), candidate_ids)
        if length(cat_idx[k]) > 0
            push!(inclusion_probs, sampler.inclusion_dict[categories[k]])
        else
            push!(inclusion_probs, 0)
        end
    end

    inclusion_probs ./= sum(inclusion_probs)

    # check if there are empty categories, if so set incl prob to 0 and
    # renormalize?
    selected_cats = rand(Categorical(inclusion_probs), n)
    for (i,c) in enumerate(selected_cats)
        if length(cat_idx[c]) > 0 
            coords[i] = candidates[rand(cat_idx[c])]
        end 
    end
    return coords
end

# Utils
_default_pool(::SpatiallyStratified) = pool(_default_strata((50,50)))
_default_categories() = Dict{Int,String}(1=>"A", 2=>"B", 3=>"C")
_default_inclusion() = Dict{String,Float64}("A"=>0.5, "B"=>0.3, "C"=>0.2)

function _default_strata(sz)
    mat = zeros(typeof(sz[1]), sz...)

    x = sz[1] ÷ 2
    y = sz[2] ÷ 3

    mat[begin:x, :] .= 1
    mat[(x + 1):end, begin:y] .= 2
    mat[(x + 1):end, (y + 1):end] .= 3
    return Layer(mat)
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "SpatiallyStratified default constructor works" begin
    @test typeof(SpatiallyStratified()) <: SpatiallyStratified
end

@testitem "SpatiallyStratified with default arguments can generate points" begin
    ss = SpatiallyStratified()
    coords = sample(ss)
    @test typeof(coords) <: Sites
end

@testitem "SpatiallyStratified throws error when number of sites is below 2" begin
    @test_throws TooFewSites SpatiallyStratified(numsites = -1)
    @test_throws TooFewSites SpatiallyStratified(numsites = 0)
    @test_throws TooFewSites SpatiallyStratified(numsites = 1)
end

@testitem "SpatiallyStratified can use custom number of points as keyword argument" begin
    NUM_POINTS = 42
    ss = SpatiallyStratified(; numsites = NUM_POINTS)
    @test ss.numsites == NUM_POINTS
end
