module Architectures

export CPU, CuGPU, MtlGPU, Architecture
export array_type, rng_type
export device

using CUDA
using Metal
using GPUArrays
using KernelAbstractions
using CUDA.CUDAKernels
using Metal.MetalKernels
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
    CuGPU <: Architecture
Run PlanktonIndividuals on one CUDA GPU node.
"""
struct CuGPU <: Architecture end

"""
    MtlGPU <: Architecture
Run PlanktonIndividuals on M-series GPU on Mac.
"""
struct MtlGPU <: Architecture end

device(::CPU) = KernelAbstractions.CPU()
device(::CuGPU) = CUDABackend()
device(::MtlGPU) = MetalBackend()

array_type(::CPU) = Array
array_type(::CuGPU) = CuArray
array_type(::MtlGPU) = MtlArray

rng_type(::CPU) = MersenneTwister()
rng_type(::CuGPU) = CURAND.default_rng()
rng_type(::MtlGPU) = GPUArrays.default_rng(MtlArray)

end