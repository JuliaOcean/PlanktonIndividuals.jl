module Abiotic

export construct_abiotic_particle, initialize_abiotic_particle!
export particle_interaction!, particle_release!
export particles_from_bcs!

using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures: device, Architecture, rng_type, array_type, unsafe_free!
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Biogeochemistry: default_bcs, getbc

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode
using PlanktonIndividuals: individuals, phytoplankton, abiotic_particle
using PlanktonIndividuals: BoundaryConditions

include("../utils.jl")
include("particle_generation.jl")
include("particle_interaction.jl")

end