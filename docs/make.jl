push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterCitations
using DocumenterVitepress
using BiodiversityObservationNetworks

bib = CitationBibliography(joinpath(@__DIR__, "BONs.bib"))

makedocs(
    sitename = "BiodiversityObservationNetworks.jl",
    authors = "Michael D. Catchen, Timothée Poisot, Kari Norman, Hana Mayall, Tom Malpas",
    modules = [BiodiversityObservationNetworks],
    format = DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl",
        devurl="dev",
    ),
    warnonly = true,
    plugins = [bib]
)


deploydocs(;
    repo = "github.com/PoisotLab/BiodiversityObservationNetworks.jl.git",
    push_preview = true,
)
