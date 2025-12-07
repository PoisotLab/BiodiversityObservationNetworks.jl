# Load the current repo version of BONs.jl
bonsjl_path = dirname(dirname(Base.current_project()))
push!(LOAD_PATH, bonsjl_path)
using BiodiversityObservationNetworks

# Load the rest of the build environment
using Documenter
using DocumenterVitepress
using DocumenterCitations
using Literate
using Markdown
using InteractiveUtils
using Dates
using PrettyTables

const bibfile = joinpath(bonsjl_path, "docs", "BONs.bib")

# Cleanup the bibliography file to make DocumenterCitations happy 
lines = readlines(bibfile)
open(bibfile, "w") do bfile
    for line in lines
        if contains(line, "journaltitle")
            println(bfile, replace(line, "journaltitle" => "journal"))
        elseif contains(line, "date")
            yrmatch = match(r"{(\d{4})", line)
            if !isnothing(yrmatch)
                println(bfile, "year = {$(yrmatch[1])},")
            end
            println(bfile, line)
        else
            println(bfile, line)
        end
    end
end
bib = CitationBibliography(
    bibfile;
    style = :authoryear,
)


makedocs(;
    sitename = "BiodiversityObservationNetworks.jl",
    format = MarkdownVitepress(;
        repo = "github.com/PoisotLab/BiodiversityObservationNetworks.jl",
    ),
    warnonly = true,
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/PoisotLab/BiodiversityObservationNetworks.jl.git",
    devbranch = "main",
    push_preview = true,
)
