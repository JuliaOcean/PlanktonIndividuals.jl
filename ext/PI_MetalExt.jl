module PI_MetalExt

using Metal
using Metal.MetalKernels
import PlanktonIndividuals.Architectures: GPU, device, array_type, rng_type, isfunctional, unsafe_free!

device(::GPU) = MetalBackend()
array_type(::GPU) = MtlArray
rng_type(::GPU) = Metal.default_rng()
isfunctional(::GPU) = Metal.functional()
unsafe_free!(m::MtlArray) = Metal.unsafe_free!(m)

end