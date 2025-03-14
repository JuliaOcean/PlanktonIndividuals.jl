module Abiotic

using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures: device, Architecture, rng_type, array_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

include("../utils.jl")

end