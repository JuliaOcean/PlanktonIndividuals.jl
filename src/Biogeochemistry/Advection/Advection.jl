module Advection

export adv_flux_x, adv_flux_y, adv_flux_z
export tracer_advection!
export tracer_diffusion!
export tracer_sinking!

using KernelAbstractions

using PlanktonIndividuals.Grids
using PlanktonIndividuals.Architectures: device, array_type, Architecture
using PlanktonIndividuals.Biogeochemistry: tracer_names
using PlanktonIndividuals.Biogeochemistry: fill_halo_Gcs!, fill_halo_tracer!

include("operaters.jl")
include("tracer_diffusion.jl")
include("DST3FL.jl")
include("multi_dim_adv.jl")
include("tracer_sinking.jl")

end

