
"""
Group
"""
@kwdef mutable struct Group
    name::Union{Symbol,String} = :none
    id::Union{<:Integer, Nothing} = nothing
end
id(g::Group) = g.id
name(g::Group) = g.name
function Base.show(io::IO, ::MIME"text/plain", g::Group) 
    tprint(
        io,
        "::Group called {yellow}$(Symbol(name(g))){/yellow}"
    )
end 


function Base.show(io::IO, ::MIME"text/plain", gs::Vector{Group}) 
    print(
        io,
        string(
            Panel(
                highlight("""
                a set of {italic}$(length(unique(gs))) unique {/italic}groups corresponding to {italic}$(length(gs)){/italic} layers 
                 
                - {blue}unique groups{/blue} => $([name(g) for g in unique(gs)])
                
                - {red}all groups{/red} => $(map(name,gs))
                """),
                title = highlight(string(typeof(gs), "[::",typeof(gs[1]), "]")),
                style = "grey dim",
                title_style = "default green bold",
                padding = (2, 2, 1, 1),
                width = 60,
            ),
        ),
    )
end 
