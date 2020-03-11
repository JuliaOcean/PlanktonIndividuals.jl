using Documenter
using DocumenterMarkdown
using PlanktonIndividuals

makedocs(
    sitename = "PlanktonIndividuals",
    format = Documenter.HTML(),
    modules = [PlanktonIndividuals],
    pages = [
    "Home" => "index.md",
    "Various" => "various.md",
            ]
)

deploydocs(;
    repo="github.com/zhenwu0728/PlanktonIndividuals.jl.git/",
)
