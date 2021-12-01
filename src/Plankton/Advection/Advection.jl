module Advection

export plankton_advection!
export plankton_diffusion!

export get_xf_index, get_yf_index, get_zf_index

using KernelAbstractions
using CUDA
using Random

using PlanktonIndividuals.Architectures: device, Architecture, rng_type
using PlanktonIndividuals.Grids

include("interpolation.jl")
include("plankton_advection_kernels.jl")
include("plankton_advection.jl")
include("plankton_diffusion.jl")

end