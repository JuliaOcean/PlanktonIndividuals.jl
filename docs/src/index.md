# PhytoAgentModel.jl

Documentation for `PhytoAgentModel.jl` which simulates the behavior of an ensemble of phytoplankton `agents`.

### Use Example

Here we use [Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl) to generate velocity fields and then use those to drive the agent-based model.

```
Pkg.develop(PackageSpec(path="AgentPhytModel_3D"))
using PhytoAgentModel
p = dirname(pathof(PhytoAgentModel))
include(joinpath(p,"Oceananigans_PlanktonAgents.jl"))
```

### Unit Testing

The tests use input files from `AgentPhytModel_3D/samples/` and then compare the results to `samples/testB1B2*.csv`. The test suite includes zero-, one-, two-, and three-dimensional simulations.

```
Pkg.develop(PackageSpec(path="AgentPhytModel_3D"))
Pkg.test("PhytoAgentModel")
```

_Contents:_

```@contents
Pages = ["index.md","various.md"]
Depth = 3
```

## API Guide

```@index
```

```@autodocs
Modules = [PhytoAgentModel]
Order   = [:type,:function]
```

