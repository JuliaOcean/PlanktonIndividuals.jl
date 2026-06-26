module ColonyIronEnergy

export colony_update!
export construct_colony, initialize_colony!

using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures: device, Architecture, rng_type, array_type, unsafe_free!
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode
using PlanktonIndividuals: individuals, phytoplankton, colony_particle, abiotic_particle

include("../../utils.jl")
include("../../Plankton/IronEnergyMode/growth_kernels.jl")
include("../../Plankton/division_death_probability.jl")
include("../../Plankton/IronEnergyMode/consume_loss.jl")
include("../../Plankton/IronEnergyMode/division_death.jl")
include("colony_generation.jl")
include("material_exchange.jl")
include("colony_growth.jl")
include("colony_update.jl")

end
