# PlanktonIndividuals.jl
[![Build Status](https://travis-ci.org/JuliaOcean/PlanktonIndividuals.jl.svg?branch=master)](https://travis-ci.org/JuliaOcean/PlanktonIndividuals.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/dev)
[![DOI](https://zenodo.org/badge/178023615.svg)](https://zenodo.org/badge/latestdoi/178023615)

This package simulates the behavior of an ensemble of phytoplankton `individuals`.

### Use Example

Here we use [Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl) to generate velocity fields and then use those to drive the individual-based model.

```
Pkg.develop(PackageSpec(path="PlanktonIndividuals.jl"))
using PlanktonIndividuals
p = dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"examples/Oceananigans_PlanktonIndividuals.jl"))
```

### Unit Testing
The tests use input files from `samples/`. The test suite includes zero-, one-, two-, and three-dimensional simulations.

```
Pkg.develop(PackageSpec(path="PlanktonIndividuals.jl"))
Pkg.test("PlanktonIndividuals")
```



