# Library

The public user interface.

## Architectures

```@autodocs
Modules = [PlanktonKernels.Architectures]
Private = false
Pages   = [ "Architectures.jl"]
```

## Grids

```@autodocs
Modules = [PlanktonKernels.Grids]
Private = false
Pages   = [
    "Grids/Grids.jl",
    "Grids/rectilinear_grid.jl",
    "Grids/lat_lon_grid.jl",
]
```


## Biogeochemistry

Tracer initialization and updates use `PlanktonKernels.Biogeochemistry` directly.
Use `bgc_tracer_init`, `generate_bgc_tracers`, and `bgc_tracer_update!` from that
module. Field utilities and halo operations are in `PlanktonKernels.Fields`.
Use `PlanktonKernels.Fields.set_bc!` for tracer boundaries and
`PlanktonIndividuals.set_bc_particle!` for particle boundaries.

## Parameter defaults and updates

```@autodocs
Modules = [PlanktonIndividuals.Model, PlanktonIndividuals.Simulation,
           PlanktonIndividuals.Individuals]
Public = true
Private = true
Pages = ["plankton_params.jl", "Model/models.jl", "default_forcing.jl"]
```

## Diagnostics

```@autodocs
Modules = [PlanktonIndividuals.Diagnostics]
Private = false
Pages   = [
    "Diagnostics/Diagnostics.jl",
    "Diagnostics/diagnostics_struct.jl"
]
```

## Model

```@autodocs
Modules = [PlanktonIndividuals, PlanktonIndividuals.Model]
Private = false
Pages   =[
    "PlanktonIndividuals.jl",
    "Model/Model.jl",
    "model_structs.jl",
    "Model/models.jl"
]
```

## Simulation

```@autodocs
Modules = [PlanktonIndividuals.Simulation]
Private = false
Pages   =[
    "Simulation/Simulation.jl",
    "Simulation/simulations.jl",
    "Simulation/update.jl",
    "Simulation/utils.jl"
]
```

## Output

```@autodocs
Modules = [PlanktonIndividuals.Output]
Private = false
Pages   =[
    "Output/Output.jl",
    "Output/output_writers.jl"
]
```
