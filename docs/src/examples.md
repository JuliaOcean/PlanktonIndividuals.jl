# [Examples](@id examples)

- [Lab Experiment](examples/0D_experiment.md)
- [Ocean Transect](examples/vertical_2D_example.md)
- [Two-Dimensional Map](examples/horizontal_2D_example.md)
- [Three-Dimensional Domain](examples/surface_mixing_3D_example.md)
- [Two-Dimensional Global Ocean](examples/global_ocean_2D_example.md)
- [Three-Dimensional Global Ocean](examples/global_ocean_3D_example.md)

## [Running the examples](@id examples-running)

The examples are Julia scripts in the `examples/` directory. From a local
checkout, use Julia 1.11 or newer and run:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/0D_experiment.jl
```

## Building the example pages

From the checkout root, build the documentation with:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Literate generates the Markdown pages and Documenter runs the examples while
building the site in `docs/build/`. The global 3D page shows a pre-rendered
animation by default. To regenerate its plots, use:

```sh
PLANKTON_DOCS_RUN_GLOBAL_3D=true julia --project=docs docs/make.jl
```

This run downloads about 550 MB of OCCA climatology data from the Harvard
Dataverse on the first run and reuses the local files later. In GitHub Actions,
start the **Documentation** workflow manually and enable **Run the global 3D
example** to use the same option and restore the OCCA cache.
