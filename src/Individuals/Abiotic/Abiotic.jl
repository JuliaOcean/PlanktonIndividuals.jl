module Abiotic


using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonKernels.Architectures: device, Architecture, rng_type, array_type, unsafe_free!
using PlanktonKernels.Grids: AbstractGrid, ΔxC, ΔyC, ΔzC, ΔxF, ΔyF, ΔzF, Ax, Ay, Az, volume
using PlanktonKernels.Fields: default_bcs, getbc, BoundaryConditions

using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals: AbioticParticle


include("../utils.jl")
include("particle_generation.jl")
include("particle_interaction.jl")

end
