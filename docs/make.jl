using Documenter, NeutralLandscapes

# For GR docs bug
ENV["GKSwstype"] = "100"

makedocs(
    sitename="BiodiversityObservationNetworks",
    authors="Timothée Poisot + others",
    modules=[BiodiversityObservationNetworks],
    pages=[
        "Index" => "index.md",
    ],
    checkdocs=:all,
    #strict=true,
)

deploydocs(
    deps=Deps.pip("pygments", "python-markdown-math"),
    repo="github.com/EcoJulia/NeutralLandscapes.jl.git",
    devbranch="main",
    push_preview=true
)
