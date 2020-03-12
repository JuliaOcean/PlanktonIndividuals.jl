using Documenter
using PlanktonIndividuals

makedocs(;
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(),
         pages = [
             "Home" => "index.md",
             "Various" => "various.md",
         ],
         repo="https://github.com/zhenwu0728/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
         assets=String[],
)

deploydocs(;
    repo="github.com/zhenwu0728/PlanktonIndividuals.jl.git/",
)
