using Documenter, Literate
using PlanktonIndividuals

examples = ["vertical_2D_example.jl", "horizontal_2D_example.jl"]

for i in 1:length(examples)
    INPUT = joinpath(@__DIR__, "..", "examples", examples[i])
    OUTPUT = joinpath(@__DIR__, "src","generated")

    Literate.markdown(INPUT, OUTPUT, documenter = true)
    # Literate.notebook(INPUT, OUTPUT, execute = true)
end

makedocs(;
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(),
         pages = [
             "Home" => "home.md",
            #  "Various" => "various.md",
             "Equations" => "equations.md",
             "generated/vertical_2D_example.md",
             "generated/horizontal_2D_example.md",
             ],
         repo="https://github.com/JuliaOcean/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
)

deploydocs(repo="github.com/JuliaOcean/PlanktonIndividuals.jl.git")
