# PlanktonIndividuals.jl

## Overview

`PlanktonIndividuals.jl` is a fast individual-based model written in Julia that can be run on both CPU and GPU. It simulates the life cycle of phytoplankton cells as Lagrangian particles in the ocean while nutrients are represented and advected as Eulerian tracers. The model is used to simulate and interpret the temporal and spatial variations in phytoplankton cell density, stoichiometry, as well as growth and division behaviors induced by diel cycle and physical motions ranging from sub-mesoscale to large scale processes.

`PlanktonIndividuals.jl` can simulate multiple functional groups of phytoplankton with different growth and division strategies which will illustrate the interactions within and between functional groups. The individuals can be simulated in 0-dimensional domain (like lab experiments) and 1-3 dimensional domains where individuals will be advected by velocities provided by various models or observations.

## Getting Help

If you are interested in using `PlanktonIndividuals.jl` or are trying to figure out how to use it, please feel free to ask us questions and get in touch!  

If you're trying to set up a model then maybe you want to check out the examples. Please feel free to [open an issue](https://github.com/JuliaOcean/PlanktonIndividuals.jl/issues)
if you have any questions, comments, suggestions, etc!
