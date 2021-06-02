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
    # "Table Of Content" => "examples.md",
    "Lab Experiment" => "generated/0D_experiment.md",
    "Ocean Transect" => "generated/vertical_2D_example.md",
    "Two-Dimensional Map" => "generated/horizontal_2D_example.md",
    "Three-Dimensional Domain" => "generated/surface_mixing_3D_example.md",
]
model_description = [
    # "Table Of Content" => "model_description.md",
    "Phytoplankton Physiology" => "phyto_equations.md",
    "Biogeochemistry" => "bgc_equations.md",
]

makedocs(;
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(collapselevel = 1, mathengine = MathJax3()),
         pages = [
             "Home" => "index.md",
             "GPU Support" => "gpu_support.md",
             "Model Description" => model_description,
             "Model Configuration" => "model_setup.md",
             "Model Simulation" => "model_run.md",
             "Examples"  => example_pages,
             "Benchmarks"  => "benchmarks.md",
             "Function index" => "function_index.md"
             ],
         repo="https://github.com/JuliaOcean/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
)

deploydocs(repo="github.com/JuliaOcean/PlanktonIndividuals.jl.git")
