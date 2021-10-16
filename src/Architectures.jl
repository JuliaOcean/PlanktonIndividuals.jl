module Architectures

export CPU, GPU, Architecture
export array_type, rng_type
export device

using CUDA
using KernelAbstractions, CUDAKernels
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
Run PlanktonIndividuals on one GPU node.
"""
struct GPU <: Architecture end

device(::CPU) = KernelAbstractions.CPU()
device(::GPU) = CUDAKernels.CUDADevice()

array_type(::CPU) = Array
if CUDA.has_cuda()
    array_type(::GPU) = CuArray
end

rng_type(::CPU) = MersenneTwister()
if CUDA.has_cuda()
    rng_type(::GPU) = CURAND.default_rng()
end

end
