using Documenter, Literate, PlanktonIndividuals, PlanktonKernels

# Use a headless GR workstation for both local and CI builds.
get!(ENV, "GKSwstype", "100")

examples = [
    "Lab Experiment" => "0D_experiment",
    "Ocean Transect" => "vertical_2D_example",
    "Two-Dimensional Map" => "horizontal_2D_example",
    "Three-Dimensional Domain" => "surface_mixing_3D_example",
    "Two-Dimensional Global Ocean" => "global_ocean_2D_example",
    "Three-Dimensional Global Ocean" => "global_ocean_3D_example",
]

examples_dir = joinpath(@__DIR__, "..", "examples")
generated_dir = joinpath(@__DIR__, "src", "examples")
mkpath(generated_dir)
run_global_3d = get(ENV, "PLANKTON_DOCS_RUN_GLOBAL_3D", "false") == "true"

for (_, name) in examples
    source = joinpath(examples_dir, name * ".jl")
    if name == "global_ocean_3D_example" && !run_global_3d
        animation = "anim_3D_global.gif"
        cp(joinpath(examples_dir, "figures", animation), joinpath(generated_dir, animation); force=true)
        preview = """

        ## Pre-rendered animation

        Standard documentation builds show the code and this animation from a prior
        run. Set `PLANKTON_DOCS_RUN_GLOBAL_3D=true` to execute the example and
        generate its plots using the OCCA climatology.

        ![Global ocean animation from a prior run]($animation)
        """
        Literate.markdown(source, generated_dir; execute=false,
                          codefence="````julia" => "````",
                          postprocess=markdown -> markdown * preview)
    else
        Literate.markdown(source, generated_dir; execute=false)
    end
end

example_pages = ["Overview" => "examples.md"]
append!(example_pages, [title => "examples/$name.md" for (title, name) in examples])

model_description = [
    "Phytoplankton Physiology" => "phyto_equations.md",
    "Biogeochemistry" => "bgc_equations.md",
]

makedocs(;
         root = @__DIR__,
         modules = [PlanktonIndividuals],
         format = Documenter.HTML(collapselevel = 1, mathengine = MathJax3(), edit_link = :commit),
         pages = [
             "Home" => "index.md",
             "Model Description" => model_description,
             "Model Configuration" => "model_setup.md",
             "Model Simulation" => "model_run.md",
             "GPU Support" => "gpu_support.md",
             "Examples" => example_pages,
             "Benchmarks" => "benchmarks.md",
             "Library" => "library.md",
             "Function index" => "function_index.md",
         ],
         warnonly = [:cross_references, :missing_docs],
         doctest = false,
         repo = "https://github.com/JuliaOcean/PlanktonIndividuals.jl/blob/{commit}{path}#L{line}",
         sitename = "PlanktonIndividuals.jl",
         authors = "Zhen Wu <zhenwu@mit.edu>",
)

deploydocs(repo="github.com/JuliaOcean/PlanktonIndividuals.jl.git", devbranch="main")
