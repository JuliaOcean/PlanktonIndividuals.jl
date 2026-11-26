# PlanktonIndividuals.jl

[![Linux](https://github.com/JuliaOcean/PlanktonIndividuals.jl/actions/workflows/linux.yml/badge.svg)](https://github.com/JuliaOcean/PlanktonIndividuals.jl/actions/workflows/linux.yml)
[![doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/stable)
[![doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaOcean.github.io/PlanktonIndividuals.jl/dev)
[![codecov](https://codecov.io/gh/JuliaOcean/PlanktonIndividuals.jl/branch/master/graph/badge.svg?token=jJL053vHAM)](https://codecov.io/gh/JuliaOcean/PlanktonIndividuals.jl)
[![DOI](https://zenodo.org/badge/178023615.svg)](https://zenodo.org/badge/latestdoi/178023615)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04207/status.svg)](https://doi.org/10.21105/joss.04207)

![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_3D_global.gif)

`PlanktonIndividuals.jl` is a fast individual-based model written in Julia that can be run on both CPU and GPU. It simulates the life cycle of phytoplankton cells as Lagrangian particles in the ocean while nutrients are represented as Eulerian, density-based tracers using a [3rd order advection scheme](https://mitgcm.readthedocs.io/en/latest/algorithm/adv-schemes.html#third-order-direct-space-time-with-flux-limiting). The model is used to simulate and interpret the temporal and spacial variations of phytoplankton cell densities and stoichiometry as well as growth and division behaviors induced by diel cycle and physical motions ranging from sub-mesoscale to large scale processes.

### Developing with PlanktonKernels

Shared architectures, grids, fields, transport, and biogeochemical numerics are
provided by PlanktonKernels. With sibling checkouts, initialize the local dependency:

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path="../PlanktonKernels.jl")
Pkg.instantiate()
Pkg.test()
```

Import shared grids, architectures, units, and biogeochemical helpers directly
from `PlanktonKernels`. Field-copy helpers such as `vel_copy!` are available from
`PlanktonKernels.Fields`.

Tracer initialization now uses `PlanktonKernels.Biogeochemistry` directly.
Use `PlanktonKernels.Fields.set_bc!` for tracer boundaries and
`PlanktonIndividuals.set_bc_particle!` for particle boundaries. GPU methods are supplied by PlanktonKernels extensions
when CUDA or Metal is loaded. Use one GPU backend per session.

Individual defaults are defined in `src/Individuals/plankton_params.jl`, with
separate methods for each physiology mode. Particle diffusivities `κhP` and `κvP`
are per-species parameter arrays (m²/s), defaulting to zero, supplied through
`phyto_setup`, `abiotic_setup`, or `colony_setup` parameters. They are no longer
biogeochemical parameters. A colony moves as one particle using the diffusivities
of its first species; other species follow that position.


Model-specific settings live in `model.options`:

```julia
opt = model_options()
opt.kc = 0.05
opt.kw = 0.046
opt.shared_graz = 0.0
opt.max_individuals = 16384
model = PlanktonModel(CPU(), grid; options=opt)
```

`bgc_params` now accepts only PlanktonKernels biogeochemical parameters; light
attenuation and grazing settings belong in options. Set particle limits through `options.max_individuals` and
`options.max_candidates`. Each model retains the supplied options object. Particle limits determine allocations
at construction and should not be changed afterward.
