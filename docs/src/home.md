# PlanktonIndividuals.jl

Documentation for `PlanktonIndividuals.jl` which simulates the behavior of an ensemble of phytoplankton `individuals`.

## Introduction

Here's the introduction of `PlanktonIndividuals.jl`

<!-- ## Use Example

Here we use [Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl) to generate velocity fields and then use those to drive the individual-based model.

```
Pkg.develop(PackageSpec(path="PlanktonIndividuals.jl"))
using PlanktonIndividuals
p = dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"Oceananigans_PlanktonIndividuals.jl"))
``` -->

## Unit Testing

The tests use input files from `samples/`. The test suite includes zero-, one-, two-, and three-dimensional simulations.

```
Pkg.develop(PackageSpec(path="PlanktonIndividuals.jl"))
Pkg.test("PlanktonIndividuals")
```

_Contents:_

```@contents
Pages = ["home.md", "equations.md"]
Depth = 3
```

## API Guide

```@index
```

```@autodocs
Modules = [PlanktonIndividuals]
Order   = [:type,:function]
```

