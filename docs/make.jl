push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterCitations
using DocumenterMarkdown
using BiodiversityObservationNetworks

bibliography = CitationBibliography(joinpath(@__DIR__, "BONs.bib"))

makedocs(
    bibliography;
    sitename = "BiodiversityObservationNetwork.jl",
    authors = "Michael D. Catchen, Timothée Poisot, Kari Norman, Hana Mayall, Tom Malpas",
    modules = [BiodiversityObservationNetworks],
    format = Markdown(),
    
)


deploydocs(;
    deps = Deps.pip("mkdocs", "pygments", "python-markdown-math", "mkdocs-material"),
    repo = "github.com/PoisotLab/BiodiversityObservationNetworks.jl.git",
    devbranch = "main",
    make = () -> run(`mkdocs build`),
    target = "site",
    push_preview = true,
)
