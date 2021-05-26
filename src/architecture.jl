abstract type Architecture end

struct CPU <: Architecture end

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
