module Architectures

export CPU, GPU, Architecture
export array_type, rng_type
export device, isfunctional
export unsafe_free!

using KernelAbstractions

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

##### CPU #####
device(::CPU) = KernelAbstractions.CPU()
array_type(::CPU) = Array
rng_type(::CPU) = MersenneTwister()
isfunctional(::CPU) = true
unsafe_free!(m::Array) = nothing

end