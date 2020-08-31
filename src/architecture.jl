abstract type Architecture end

struct CPUs <: Architecture end

struct GPUs <: Architecture end

macro hascuda(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end

device(::CPUs) = KernelAbstractions.CPU()
device(::GPUs) = KernelAbstractions.CUDADevice()

array_type(::CPUs) = Array
@hascuda array_type(::GPUs) = CuArray
