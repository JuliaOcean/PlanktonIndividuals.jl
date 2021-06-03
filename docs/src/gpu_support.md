# GPU Support

`PlanktonIndividuals.jl` has support from `CUDA.jl` and `KernelAbstractions.jl` to be able to run on graphical processing unit (GPU) for higher performance. Depending on the combination of CPU and GPU you have, a speedup of 35x is possible. Please see [Benchmarks](@ref benchmarks) for more details.

## How to use a GPU

To use a GPU to run `PlanktonIndividuals.jl` is easy. Users do not need to rewrite the setup or simulation script to change the architecture to run on. See [Architecture](@ref) for detailed instructions on setting up a model on GPU.

!!! tip "Running on GPUs"
    If you are having issues with running `PlanktonIndividuals` on a GPU, please
    [open an issue](https://github.com/JuliaOcean/PlanktonIndividuals.jl/issues/new)

## When to use a GPU

GPU is very useful when running large simulations (either large domain or huge number of individuals, or both). If you simulate over 10,000 individuals, you will probably benefit form GPU. Please note, GPU is usually memory-limited, that is to say, you will probably fill up the memory on GPU long before the model slows down.

`Individuals` take up a large amount of GPU memory due to complicated physiological processes and diagnostic demand. For now, please do not try more than 50,000 individuals for a 12GB GPU.

## GPU resources

There are a few resources you can try to acquire a GPU from.

1. Google Colab provides GPUs but you need to install Julia manually. Please see [this post](https://discourse.julialang.org/t/julia-on-google-colab-free-gpu-accelerated-shareable-notebooks/15319/39) on the Julia Discourse for detailed instructions.

2. [Code Ocean](https://codeocean.com/) also has [GPU support](https://help.codeocean.com/en/articles/1053107-gpu-support). You can use "Ubuntu Linux with GPU support (18.04.3)" but you still have to install Julia manually.
