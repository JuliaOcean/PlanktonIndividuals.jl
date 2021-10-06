module Plankton

export plankton_advection!
export plankton_diffusion!
export plankton_update!
export individuals, gen_individuals!
export find_inds!, find_NPT!, acc_counts!, calc_par!

using CUDA
using StructArrays
using Random
using KernelAbstractions: @index, @kernel, Event, MultiEvent, wait
using KernelAbstractions.Extras.LoopInfo: @unroll

using PlanktonIndividuals.Architectures: device, Architecture, GPU, CPU, rng_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Output

include("Advection/Advection.jl")
include("Quota/Quota.jl")
include("utils.jl")

using .Advection
using .Quota

end