![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_global.gif)

# PlanktonIndividuals.jl

[![Build Status](https://travis-ci.org/JuliaOcean/PlanktonIndividuals.jl.svg?branch=master)](https://travis-ci.org/JuliaOcean/PlanktonIndividuals.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/dev)
[![codecov](https://codecov.io/gh/JuliaOcean/PlanktonIndividuals.jl/branch/master/graph/badge.svg?token=jJL053vHAM)](https://codecov.io/gh/JuliaOcean/PlanktonIndividuals.jl)
[![DOI](https://zenodo.org/badge/178023615.svg)](https://zenodo.org/badge/latestdoi/178023615)

This package simulates the behavior of an ensemble of phytoplankton _individuals_.

## Use Examples

### 1. Simple Flow Fields In Two Dimensions

```julia
using PlanktonIndividuals
p = dirname(pathof(PlanktonIndividuals))
#include(joinpath(p,"../examples/vertical_2D_example.jl"))
include(joinpath(p,"../examples/horizontal_2D_example.jl"))
```

### 2. Closer Look Into One Grid Box

```julia
using PlanktonIndividuals
p = dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"../examples/0D_experiment.jl"))
```

### 3. Turbulent Flow Fields In Three Dimensions

Here [Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl) is used to generate velocity fields and then use those to drive the individual-based model.

```julia
using PlanktonIndividuals
p = dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"../examples/surface_mixing_3D_example.jl"))
```

## Unit Testing

The test suite includes zero-, one-, two-, and three-dimensional simulations using input files from `samples/`.

```julia
using Pkg; Pkg.test("PlanktonIndividuals")
```

## Installation

To add `PlanktonIndividuals.jl` to your Julia environment:

```julia
using Pkg; Pkg.add("PlanktonIndividuals.jl")
```

or if you cloned the repository via `git` directly:

```julia
using Pkg; 
p=PackageSpec(path="PlanktonIndividuals.jl")
Pkg.develop()
```
