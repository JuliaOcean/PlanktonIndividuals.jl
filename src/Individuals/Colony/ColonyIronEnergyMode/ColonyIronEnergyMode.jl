module ColonyIronEnergy


using KernelAbstractions
using StructArrays
using Random
using LinearAlgebra: dot

using PlanktonKernels.Architectures: device, Architecture, rng_type, array_type, unsafe_free!
using PlanktonKernels.Grids: AbstractGrid, ΔzF, volume

using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Diagnostics: diags_proc!, diags_colony!
using PlanktonIndividuals: IronEnergyMode, Phytoplankton, ColonyParticle


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
