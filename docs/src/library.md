# Library

The public user interface.

## Architectures

```@autodocs
Modules = [PlanktonKernels.Architectures]
Private = false
Pages   = [ "Architectures.jl"]
```

## Grids

See [Model Configuration](model_setup.md#Grid) for regular and latitude–longitude grids.

## Biogeochemistry

Tracer initialization and updates use `PlanktonKernels.Biogeochemistry` directly.
Use `bgc_tracer_default_init`, `generate_bgc_tracers`, and `bgc_tracer_update!` from that
module. Field utilities and halo operations are in `PlanktonKernels.Fields`.
Use `PlanktonKernels.Fields.set_bc!` for tracer boundaries and
`PlanktonIndividuals.set_bc_particle!` for particle boundaries.

## Parameter defaults and updates

```@autodocs
Modules = [PlanktonIndividuals.Simulation, PlanktonIndividuals.Individuals]
Public = true
Private = true
Pages = ["plankton_params.jl", "Simulation/simulations.jl"]
Filter = object -> object in (default_PARF, default_temperature, phyt_params_default, colony_params_default, abiotic_params_default, update_phyt_params, update_colony_params, update_abiotic_params)
```

## Diagnostics

```@autodocs
Modules = [PlanktonIndividuals.Diagnostics]
Private = false
Filter = object -> object !== PlanktonDiagnostics
Pages   = [
    "Diagnostics/Diagnostics.jl",
    "Diagnostics/diagnostics_struct.jl"
]
```

## Model

```@autodocs
Modules = [PlanktonIndividuals, PlanktonIndividuals.Model]
Private = false
Filter = object -> object !== PlanktonModel
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
Filter = object -> !(object in (PlanktonSimulation, sim_run!, default_PARF, default_temperature))
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
Filter = object -> object !== PlanktonOutputWriter
Pages   =[
    "Output/Output.jl",
    "Output/output_writers.jl"
]
```
