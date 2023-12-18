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
    authors = "M.D. Catchen, Timothée Poisot, Kari Norman, Hana Mayall, Tom Malpas",
    modules = [BiodiversityObservationNetworks],
    pages = [
        "Index" => "index.md",
        "Vignettes" => [
            "Overview" => "vignettes/overview.md",
            "Entropy" => "vignettes/entropize.md",
            "Environmental uniqueness" => "vignettes/uniqueness.md",
        ],
    ],
    checkdocs = :all,
    warnonly = true,
    format = HTML(;size_threshold= nothing)
)

deploydocs(;
    deps = Deps.pip("pygments", "python-markdown-math"),
    repo = "github.com/EcoJulia/BiodiversityObservationNetworks.jl.git",
    devbranch = "main",
    push_preview = true,
)
