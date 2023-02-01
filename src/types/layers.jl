

"""
Layer

"""
@kwdef struct Layer
    name::Union{Symbol,String}
    layer = rand(50,100)
    group::Group
end
function Base.show(io::IO, ::MIME"text/plain", l::Layer) 
tprint(
    io,
    "::Layer called {yellow}$(Symbol(name(l))){/yellow} that is $(size(l.layer)[1]) px tall and $(size(l.layer)[2]) px wide"
)
end 
name(l::Layer) = l.name
getgroup(l::Layer) = l.group
getlayer(l::Layer) = l.layer
Base.size(l::Layer) = size(getlayer(l))
Base.keys(l::Layer) = eachindex(getlayer(l))


 # LayerSet 
 @kwdef struct LayerSet
    layers::Vector{Layer} = [Layer() for _ in 1:5]
    targets::Vector{Target} = [Target() for _ in 1:5]
end 
function Base.show(io::IO, ::MIME"text/plain", ls::LayerSet) 
    print(
        io,
        string(
            Panel(
                highlight("""
                set of {italic}$(numlayers(ls)) layers{/italic} with {italic}$(numtargets(ls)) targets{/italic} and {italic}$(numgroups(ls)) layer{/italic} groups.


                 - 🗺️  {green}layers{/green} => $([name(l) for l in getlayers(ls)])
                 
                 - 📌 {red}targets{/red} => $([name(t) for t in gettargets(ls)])

                 - 🗂️  {blue}groups{/blue} => $([name(g) for g in unique(getgroups(ls))])
                
                """),
                title = ("{yellow}::LayerSet{/yellow}"),
                style = "grey dim",
                title_style = "default green bold",
                padding = (2, 2, 1, 1),
                width = 60,
            ),
        ),
    )
end 

getlayers(ls::LayerSet) = ls.layers
getlayer(ls::LayerSet, i::Integer) = getlayers(ls)[i]
getgroups(ls::LayerSet) = [getgroup(l) for l in getlayers(ls)]
gettargets(ls::LayerSet) = ls.targets

Base.length(ls::LayerSet) = length(getlayers(ls))
Base.iterate(ls::LayerSet, i=1) = i <= length(ls) ? (getlayers(ls)[i],i+1) : nothing

numlayers(ls::LayerSet) = length(ls)
numgroups(ls::LayerSet) = length(unique(getgroups(ls)))
numtargets(ls::LayerSet) = length(gettargets(ls))
