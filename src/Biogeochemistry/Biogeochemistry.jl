module Biogeochemistry

export generate_tracers
export tracers_init, default_tracer_init
export tracer_update!

using KernelAbstractions

using PlanktonIndividuals.Architectures: device, array_type, Architecture
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields

include("Advection/Advection.jl")

using .Advection

include("tracer_fields.jl")
include("tracer_forcings.jl")
include("tracer_update.jl")

end