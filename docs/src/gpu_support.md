# GPU Support

`PlanktonIndividuals.jl` has support from `CUDA.jl` and `KernelAbstractions.jl` to be able to run on graphical processing unit (GPU) for higher performance. Depending on the combination of CPU and GPU you have, a speedup of 35x is possible. Please see [Benchmarks](@ref benchmarks) for more details.

## How to use a GPU

To use a GPU to run `PlanktonIndividuals.jl` is easy. Users do not need to rewrite the setup or simulation script to change the architecture to run on. See [Architecture](@ref) for detailed instructions on setting up a model on GPU.

!!! tip "Running on GPUs"
    If you are having issues with running `PlanktonIndividuals` on a GPU, please
    [open an issue](https://github.com/JuliaOcean/PlanktonIndividuals.jl/issues/new)

## When to use a GPU

## GPU resources
