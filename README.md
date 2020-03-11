# PlanktonIndividuals.jl

[![Travis Build Status](https://travis-ci.org/zhenwu0728/PlanktonIndividuals.svg?branch=master)](https://travis-ci.org/zhenwu0728/PlanktonIndividuals)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://zhenwu0728.github.io/PlanktonIndividuals/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://zhenwu0728.github.io/PlanktonIndividuals/dev)

This package simulates the behavior of an ensemble of phytoplankton `individuals`.

### Use Example

Here we use [Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl) to generate velocity fields and then use those to drive the individual-based model.

```
Pkg.develop(PackageSpec(path="PlanktonIndividuals.jl"))
using PlanktonIndividuals
p = dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"Oceananigans_PlanktonIndividuals.jl"))
```

### Unit Testing
The tests use input files from `samples/`. The test suite includes zero-, one-, two-, and three-dimensional simulations.

```
Pkg.develop(PackageSpec(path="PlanktonIndividuals.jl"))
Pkg.test("PlanktonIndividuals")
```



