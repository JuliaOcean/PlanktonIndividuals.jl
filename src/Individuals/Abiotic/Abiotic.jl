module Abiotic

#export abiotic_particle_update!
export construct_abiotic_particle, initialize_abiotic_particle!
export abiotic_particle_update!

using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures: device, Architecture, rng_type, array_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode
using PlanktonIndividuals: individuals, phytoplankton, abiotic_particle

include("../utils.jl")
include("particle_generation.jl")
include("particle_update.jl")

end