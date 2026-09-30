module IronEnergy


using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonKernels.Architectures: device, Architecture, rng_type, array_type, unsafe_free!
using PlanktonKernels.Grids: AbstractGrid, ΔzF, volume

using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Diagnostics: diags_spcs!, diags_proc!
using PlanktonIndividuals: Phytoplankton


include("../../utils.jl")
include("../division_death_probability.jl")
include("plankton_generation.jl")
include("growth_kernels.jl")
include("plankton_growth.jl")
include("consume_loss.jl")
include("division_death.jl")
include("plankton_update.jl")

end
