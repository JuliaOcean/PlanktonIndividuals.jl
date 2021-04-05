abstract type Architecture end

struct CPU <: Architecture end

struct GPU <: Architecture end

macro hascuda(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end

device(::CPU) = KernelAbstractions.CPU()
device(::GPU) = KernelAbstractions.CUDADevice()

array_type(::CPU) = Array
@hascuda array_type(::GPU) = CuArray

rng_type(::CPU) = MersenneTwister()
@hascuda rng_type(::GPU) = CURAND.default_rng()
