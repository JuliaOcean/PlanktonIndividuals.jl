using Documenter, Literate
using PlanktonIndividuals

EXAMPLE = joinpath(@__DIR__, "..", "examples", "vertical_2D_example.jl")
OUTPUT = joinpath(@__DIR__, "src","generated")

Literate.markdown(EXAMPLE, OUTPUT, documenter = true)
Literate.notebook(EXAMPLE, OUTPUT, execute = true)

makedocs(;
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(),
         pages = [
             "Home" => "index.md",
             "Various" => "various.md",
             "Equations" => "equations.md",
             "generated/vertical_2D_example.md",
             ],
         repo="https://github.com/JuliaOcean/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
)

deploydocs(repo="github.com/JuliaOcean/PlanktonIndividuals.jl.git")
