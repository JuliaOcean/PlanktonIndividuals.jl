using Documenter, Literate
using PlanktonIndividuals

examples = ["vertical_2D_example.jl", "horizontal_2D_example.jl", "surface_mixing_3D_example.jl", "0D_experiment.jl"]

for i in 1:length(examples)
    INPUT = joinpath(@__DIR__, "..", "examples", examples[i])
    OUTPUT = joinpath(@__DIR__, "src","generated")

    Literate.markdown(INPUT, OUTPUT, documenter = true)
    # Literate.notebook(INPUT, OUTPUT, execute = true)
end

example_pages = [
    "Vertical Two-Dimensional Example" => "generated/vertical_2D_example.md",
    "Horizontal Two-Dimensional Example" => "generated/horizontal_2D_example.md",
    "Surface Mixing Three-Dimensional Example" => "generated/surface_mixing_3D_example.md",
    "Zero-Dimensional Experiment" => "generated/0D_experiment.md"
]

makedocs(;
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(),
         pages = [
             "Home" => "index.md",
            #  "Various" => "various.md",
             "Equations" => "equations.md",
             "Examples"  => example_pages,
             "Benchmarks"  => "benchmarks.md",
             "Function index" => "function_index.md"
             ],
         repo="https://github.com/JuliaOcean/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
)

deploydocs(repo="github.com/JuliaOcean/PlanktonIndividuals.jl.git")
