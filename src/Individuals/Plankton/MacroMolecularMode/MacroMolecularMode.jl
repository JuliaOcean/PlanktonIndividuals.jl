module MacroMolecular

export plankton_update!
export construct_plankton, initialize_plankton!

using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures: device, Architecture, rng_type, array_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

include("../../utils.jl")
include("../division_death_probability.jl")
include("plankton_generation.jl")
include("growth_kernels.jl")
include("plankton_growth.jl")
include("consume_loss.jl")
include("division_death.jl")
include("plankton_update.jl")

end
