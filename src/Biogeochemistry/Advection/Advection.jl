module Advection

export adv_flux_x, adv_flux_y, adv_flux_z
export tracer_advection!
export tracer_diffusion!

using KernelAbstractions

using PlanktonIndividuals.Grids
using PlanktonIndividuals.Architectures: device, array_type, Architecture
using PlanktonIndividuals.Fields

include("operaters.jl")
include("tracer_diffusion.jl")
include("DST3FL.jl")
include("multi_dim_adv.jl")

end

