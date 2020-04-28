using Documenter
using PlanktonIndividuals

makedocs(;
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(),
         pages = [
             "Home" => "index.md",
             "Various" => "various.md",
             "Equations" => "equations.md",
         ],
         repo="https://github.com/zhenwu0728/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
)

deploydocs(repo="github.com/zhenwu0728/PlanktonIndividuals.jl.git")
