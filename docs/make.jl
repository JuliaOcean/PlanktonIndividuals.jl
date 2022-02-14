using Documenter, Pkg, PlanktonIndividuals, Plots
import PlutoSliderServer
Pkg.precompile()

examples = ["vertical_2D_example.jl", "horizontal_2D_example.jl", "surface_mixing_3D_example.jl",
            "0D_experiment.jl", "global_ocean_2D_example.jl", "global_ocean_3D_example.jl"]

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
             "Model Description" => model_description,
             "Model Configuration" => "model_setup.md",
             "Model Simulation" => "model_run.md",
             "GPU Support" => "gpu_support.md",
             "Examples"  => "examples.md",
             "Benchmarks"  => "benchmarks.md",
             "Library" => "library.md",
             "Function index" => "function_index.md"
             ],
         repo="https://github.com/JuliaOcean/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors="ZhenWu <zhenwu@mit.edu>",
)

for i in examples
    fil_in=joinpath(@__DIR__,"..", "examples",i)
    fil_out=joinpath(@__DIR__,"build", "examples",i[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    cp(fil_in,fil_out[1:end-4]*"jl")
end

deploydocs(repo="github.com/JuliaOcean/PlanktonIndividuals.jl.git")
