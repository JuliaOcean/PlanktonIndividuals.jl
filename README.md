# PlanktonIndividuals.jl

[![Build Status](https://travis-ci.org/JuliaOcean/PlanktonIndividuals.jl.svg?branch=master)](https://travis-ci.org/JuliaOcean/PlanktonIndividuals.jl)
[![doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/stable)
[![doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/dev)
[![codecov](https://codecov.io/gh/JuliaOcean/PlanktonIndividuals.jl/branch/master/graph/badge.svg?token=jJL053vHAM)](https://codecov.io/gh/JuliaOcean/PlanktonIndividuals.jl)
[![DOI](https://zenodo.org/badge/178023615.svg)](https://zenodo.org/badge/latestdoi/178023615)

![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_global.gif)

`PlanktonIndividuals.jl` is a fast individual-based model written in Julia that can be run on both CPU and GPU. It simulates the life cycle of phytoplankton cells as Lagrangian particles in the ocean while nutrients are represented as Eulerian, density-based tracers using a [3rd order advection scheme](https://mitgcm.readthedocs.io/en/latest/algorithm/adv-schemes.html#third-order-direct-space-time-with-flux-limiting). The model is used to simulate and interpret the temporal and spacial variations of phytoplankton cell densities and stoichiometry as well as growth and division behaviors induced by diel cycle and physical motions ranging from sub-mesoscale to large scale processes.

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
