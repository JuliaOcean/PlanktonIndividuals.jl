module Biogeochemistry

export generate_tracers
export tracers_init, default_tracer_init
export tracer_update!
export Field
export interior, zero_fields!
export fill_halo_vel!
export default_bcs, getbc, apply_bcs!
export set_bc!, set_bc_particle!, validate_bcs
export tracer_names

using KernelAbstractions

using PlanktonIndividuals.Architectures: device, array_type, Architecture
using PlanktonIndividuals.Grids

using PlanktonIndividuals: BoundaryConditions

const tracer_names=(:DIC,:NH4,:NO3,:PO4,:FeT,:DOC,:DON,:DOP,:DOFe,:POC,:PON,:POP,:POFe)

include("HaloRegions/halo_regions.jl")
include("HaloRegions/boundary_conditions.jl")
include("HaloRegions/apply_bcs.jl")

include("Advection/Advection.jl")

using .Advection

include("tracer_fields.jl")
include("tracer_forcings.jl")
include("tracer_update.jl")

end