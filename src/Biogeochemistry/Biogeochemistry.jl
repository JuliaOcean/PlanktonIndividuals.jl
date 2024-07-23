module Biogeochemistry

export generate_nutrients
export nutrients_init, default_nut_init
export nut_update!

using KernelAbstractions

using PlanktonIndividuals.Architectures: device, array_type, Architecture
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields

include("Advection/Advection.jl")

using .Advection

include("nutrient_fields.jl")
include("nutrient_forcings.jl")
include("nutrient_update.jl")

end