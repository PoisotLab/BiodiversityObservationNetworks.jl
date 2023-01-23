
@kwdef mutable struct Target
    name::Union{Symbol,String} = :none
    id::Union{<:Integer,Nothing} = nothing
end 
id(t::Target) = t.id
name(t::Target) = t.name
function Base.show(io::IO, ::MIME"text/plain", t::Target) 
    tprint(
        io,
        "::Target called {yellow}$(Symbol(name(t))){/yellow}"
    )
end 

function Base.show(io::IO, ::MIME"text/plain", ts::Vector{Target}) 
    tprint("a set of {italic}$(length(ts)){/italic} targets $(map(name,ts))")
end 