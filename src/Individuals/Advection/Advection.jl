module Advection

export particle_advection!
export particle_diffusion!

export get_xf_index, get_yf_index, get_zf_index

using KernelAbstractions
using Random

using PlanktonIndividuals.Architectures: device, Architecture, rng_type
using PlanktonIndividuals.Grids

include("interpolation.jl")
include("particle_advection_kernels.jl")
include("particle_advection.jl")
include("particle_diffusion.jl")

end