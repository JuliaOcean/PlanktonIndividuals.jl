# Library

The public user interface.

## Architectures

```@autodocs
Modules = [PlanktonIndividuals.Architectures]
Private = false
Pages   = [ "Architectures.jl"]
```

## Grids

```@autodocs
Modules = [PlanktonIndividuals.Grids]
Private = false
Pages   = [
    "Grids/Grids.jl",
    "Grids/regular_rectilinear_grid.jl",
    "Grids/regular_lat_lon_grid.jl",
    "Grids/vertically_stretched_lat_lon_grid.jl"
]
```

## Fields

```@autodocs
Modules = [PlanktonIndividuals.Fields]
Private = false
Pages   =[
    "Fields/Fields.jl",
    "Fields/boundary_conditions.jl"
]
```

## Biogeochmeistry

```@autodocs
Modules = [PlanktonIndividuals.Biogeochemistry]
Private = false
Pages   = [
    "Biogeochemistry/Biogeochemistry.jl",
    "Biogeochemistry/nutrient_fields.jl"
]
```

## Parameters

```@autodocs
Modules = [PlanktonIndividuals.Parameters]
Private = false
Pages   = [
    "Parameters/Parameters.jl",
    "Parameters/param_default.jl",
    "Parameters/param_update.jl"
]
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
    "Output/output_writer.jl"
]
```
