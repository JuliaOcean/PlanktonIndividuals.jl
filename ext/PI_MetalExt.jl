module PI_MetalExt

using Metal
using GPUArrays
using Metal.MetalKernels
import PlanktonIndividuals.Architectures: GPU, device, array_type, rng_type, isfunctional

device(::GPU) = MetalBackend()
array_type(::GPU) = MtlArray
rng_type(::GPU) = GPUArrays.default_rng(MtlArray)
isfunctional(::GPU) = Metal.functional()

end