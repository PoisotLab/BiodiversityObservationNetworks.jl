
struct Weights
    layers::LayerSet
    weights::Matrix  
    group_mixing::Vector
    target_mixing::Vector
end
getlayers(w::Weights) = w.layers

function Base.show(io::IO, ::MIME"text/plain", weights::Weights) 
    layer_names = map(name, getlayers(weights))
    target_names = map(name, gettargets(getlayers(weights)))
    group_names = map(name ,getgroups(getlayers(weights)))


    numtargs, numgroups = numtargets(getlayers(weights)), numtargets(getlayers(weights))


    # append title lines

    # for i in number of layers, enumerate a row     
    


    str = 
    md"""
    | {blue}Group{/blue} | {blue}Layer{/blue}  | $(target_names[1])| $(target_names[2])|
    | :----------|:----------        | ---------- |:------------:|
    | $(group_names[1]) | $(layer_names[1]) | $(weights.weights[1,1]) | $(weights.weights[1,2]) |
    | $(group_names[2]) | $(layer_names[2]) | $(weights.weights[2,1]) | $(weights.weights[2,2]) |
    | $(group_names[3]) | $(layer_names[3]) | $(weights.weights[3,1]) | $(weights.weights[3,2]) |
    | $(group_names[4]) | $(layer_names[4]) | $(weights.weights[4,1]) | $(weights.weights[4,2]) |
    | $(group_names[5]) | $(layer_names[5]) | $(weights.weights[5,1]) | $(weights.weights[5,2]) |
    """
    

    tprint(str)
end 



function _allocate_weights_matrix(layers_per_group, ntargs, ngroups)

    mat = zeros(sum(layers_per_group), ntargs)
     
    for g in 1:ngroups, t in 1:ntargs
        Ig_start = sum(layers_per_group[1:g-1])+1
        Ig_end = sum(layers_per_group[1:g])
        mat[Ig_start:Ig_end,t] .= rand(Dirichlet(ones(layers_per_group[g])))
    end
    mat
end 


function Weights(
    layers;
    target_mixing_weights = [1/numtargets(layers) for _ in 1:numtargets(layers)],
    group_mixing_weights = [1/numgroups(layers) for _ in 1:numgroups(layers)],
)
    ngroups, ntargs = numgroups(layers), numtargets(layers)
    mat = _allocate_weights_matrix(_layers_per_group(layers), ntargs, ngroups)
    Weights(layers,mat,target_mixing_weights,group_mixing_weights)
end 

_layers_per_group(layers) = begin
    uniquegroups = unique(getgroups(layers))
    groupids = [findfirst(x->getgroup(l)==x, uniquegroups) for (i,l) in enumerate(layers)]
    return StatsBase.counts(groupids)
end

function Weights(
        layers,
        mat, 
        target_mixing_weights = [1/numtargets(layers) for _ in 1:numtargets(layers)],
        group_mixing_weights = [1/numgroups(layers) for _ in 1:numgroups(layers)],
    )
    # construct adjusted mat 
    Weights(
        layers, 
        mat, 
        group_mixing_weights,
        target_mixing_weights
    )
end

rmat() = rand(50,50)

#=
layers = [Layer(rmat(), :l1), Layer(rmat(), :l2), Layer(rmat(), :l2), Layer(rmat(), :l3)]
groups = [Group(:g1), Group(:g2)]
targs =  [Target(:t1), Target(:t2)]
Weights(layers, groups, targs) =#
