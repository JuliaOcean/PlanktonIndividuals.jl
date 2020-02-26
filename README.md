# PhytoAgentModel.jl

[![Travis Build Status](https://travis-ci.org/zhenwu0728/AgentPhytModel_3D.svg?branch=master)](https://travis-ci.org/zhenwu0728/AgentPhytModel_3D)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://zhenwu0728.github.io/AgentPhytModel_3D/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://zhenwu0728.github.io/AgentPhytModel_3D/dev)

This package simulates the behavior of an ensemble of phytoplankton `agents`.

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



