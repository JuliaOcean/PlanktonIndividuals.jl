module MacroMolecular

export plankton_update!
export construct_plankton, generate_plankton!

using KernelAbstractions: @kernel, @index, Event, MultiEvent, wait
using KernelAbstractions.Extras.LoopInfo: @unroll
using CUDA
using StructArrays
using Random
using LinearAlgebra

using PlanktonIndividuals.Architectures: device, Architecture, GPU, CPU, rng_type, array_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

include("../utils.jl")
include("../division_death_probability.jl")
include("plankton_generation.jl")
include("growth_kernels.jl")
include("plankton_growth.jl")
include("consume_loss.jl")
include("division_death.jl")
include("plankton_update.jl")

end