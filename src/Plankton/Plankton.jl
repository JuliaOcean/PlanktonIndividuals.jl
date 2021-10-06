module Plankton

export plankton_advection!
export plankton_diffusion!
export plankton_update!
export individuals, gen_individuals!
export find_inds!, find_NPT!, acc_counts!, calc_par!
export CarbonMode, QuotaMode, MacroMolecularMode, AbstractMode

using CUDA
using StructArrays
using Random
using KernelAbstractions: @index, @kernel, Event, MultiEvent, wait
using KernelAbstractions.Extras.LoopInfo: @unroll

using PlanktonIndividuals.Architectures: device, Architecture, GPU, CPU, rng_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Output

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

include("Advection/Advection.jl")
include("QuotaMode/QuotaMode.jl")
include("CarbonMode/CarbonMode.jl")
include("utils.jl")

using .Advection
using .Quota
using .Carbon

end