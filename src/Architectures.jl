module Architectures

export CPU, GPU, Architecture
export array_type, rng_type
export device

using CUDA
using GPUArrays
using KernelAbstractions
using CUDA.CUDAKernels
using Random

"""
    Architecture
Abstract type for architectures supported by PlanktonIndividuals.
"""
abstract type Architecture end

"""
    CPU <: Architecture
Run PlanktonIndividuals on one CPU node.
"""
struct CPU <: Architecture end

"""
    GPU <: Architecture
Run PlanktonIndividuals on one CUDA GPU node.
"""
struct GPU <: Architecture end


device(::CPU) = KernelAbstractions.CPU()
device(::GPU) = CUDABackend()

array_type(::CPU) = Array
array_type(::GPU) = CuArray

rng_type(::CPU) = MersenneTwister()
rng_type(::GPU) = CURAND.default_rng()

end