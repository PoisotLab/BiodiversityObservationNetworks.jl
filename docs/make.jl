push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterCitations
using DocumenterVitepress
using BiodiversityObservationNetworks

bib = CitationBibliography(joinpath(@__DIR__, "BONs.bib"), style=:authoryear)

makedocs(
    sitename = "BiodiversityObservationNetworks.jl",
    authors = "Michael D. Catchen, Timothée Poisot, Kari Norman, Hana Mayall, Tom Malpas",
    modules = [BiodiversityObservationNetworks],
    format = DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl",
        devurl="dev",
    ),
    pages = [
        "Overview" => "index.md",
        "Manual" => "manual.md",
        "Tutorials" => [],
        "How To" => [],
        "Samplers" => [
            joinpath("reference","samplers", "simplerandom.md"),
            joinpath("reference", "samplers", "balancedacceptance.md"),
            joinpath("reference", "samplers", "grts.md"),
            joinpath("reference","samplers", "cube.md"),
            joinpath("reference", "samplers", "adaptivehotspot.md"),
        ],
        "Utilities" => [
            joinpath("reference", "utilities", "spatialbalance.md"),
        ],
        "API Reference" => joinpath("reference", "api.md"),
        "Bibliography" => "bibliography.md"
    ],
    warnonly = true,
    plugins = [bib]
)


deploydocs(;
    repo = "https://github.com/PoisotLab/BiodiversityObservationNetworks.jl.git",
    push_preview = true,
)
