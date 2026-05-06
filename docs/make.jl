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

@info pwd()
@info readdir()

bibfile = joinpath("docs", "BONs.bib")

# Cleanup the bibliography file to make DocumenterCitations happy despite their
# refusal to acknowledge modern field names. The people will partu like it's
# 1971 and they will like it.
lines = readlines(bibfile)
open(bibfile, "w") do bfile
    for line in lines
        if Base.contains(line, "journaltitle")
            println(bfile, replace(line, "journaltitle" => "journal"))
        elseif Base.contains(line, "date")
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
    modules=[BiodiversityObservationNetworks],
    authors="Michael D. Catchen",
    repo="github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/{commit}{path}#{line}",
    sitename="BiodiversityObservationNetworks.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/PoisotLab/BiodiversityObservationNetworks.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    warnonly=true,
    plugins = [bib],
)

# Vitepress requires their own method for deploydocs since v0.2
DocumenterVitepress.deploydocs(;
    repo = "github.com/PoisotLab/BiodiversityObservationNetworks.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main", 
    push_preview = true,
)


