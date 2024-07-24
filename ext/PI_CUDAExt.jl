module PI_CUDAExt

using CUDA
using CUDA.CUDAKernels
import PlanktonIndividuals.Architectures: GPU, device, array_type, rng_type, isfunctional

device(::GPU) = CUDABackend()
array_type(::GPU) = CuArray
rng_type(::GPU) = CURAND.default_rng()
isfunctional(::GPU) = CUDA.functional()

end