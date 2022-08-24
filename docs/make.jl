using Documenter, BiodiversityObservationNetworks
import Literate

# For GR docs bug
ENV["GKSwstype"] = "100"

vignettes = filter(
    endswith(".jl"),
    readdir(joinpath(@__DIR__, "src", "vignettes"); join = true, sort = true),
)
for vignette in vignettes
    Literate.markdown(
        vignette,
        joinpath(@__DIR__, "src", "vignettes");
        config = Dict("credit" => false, "execute" => true),
    )
end

makedocs(;
    sitename = "BiodiversityObservationNetworks",
    authors = "Timothée Poisot + others",
    modules = [BiodiversityObservationNetworks],
    pages = [
        "Index" => "index.md",
        "Vignettes" => [
            "Overview" => "vignettes/overview.md",
        ],
    ],
    checkdocs = :all,
    #strict=true,
)

deploydocs(;
    deps = Deps.pip("pygments", "python-markdown-math"),
    repo = "github.com/EcoJulia/BiodiversityObservationNetworks.jl.git",
    devbranch = "main",
    push_preview = true,
)
