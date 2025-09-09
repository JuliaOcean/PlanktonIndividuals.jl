module PI_CUDAExt

using CUDA
using CUDA.CUDAKernels
import PlanktonIndividuals.Architectures: GPU, device, array_type, rng_type, isfunctional,  unsafe_free!

device(::GPU) = CUDABackend()
array_type(::GPU) = CuArray
rng_type(::GPU) = CURAND.default_rng()
isfunctional(::GPU) = CUDA.functional()
unsafe_free!(m::CuArray) = CUDA.unsafe_free!(m)

end