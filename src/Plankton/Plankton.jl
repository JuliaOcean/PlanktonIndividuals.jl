module Plankton

export plankton_advection!
export plankton_diffusion!
export plankton_update!
export generate_individuals, individuals
export find_inds!, find_NPT!, acc_counts!, calc_par!

using CUDA
using StructArrays
using Random
using KernelAbstractions: @index, @kernel, Event, MultiEvent, wait
using KernelAbstractions.Extras.LoopInfo: @unroll

using PlanktonIndividuals.Architectures: device, Architecture, GPU, CPU, rng_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

#####
##### generate individuals of multiple species
#####
function generate_individuals(params::Dict, arch::Architecture, Nsp, N, maxN, g::AbstractGrid, mode::AbstractMode)
    plank_names = Symbol[]
    plank_data=[]
    for i in 1:Nsp
        name = Symbol("sp"*string(i))
        plank = construct_plankton(arch, i, params, maxN, mode::AbstractMode)
        generate_plankton!(plank, N, g, arch, mode)
        push!(plank_names, name)
        push!(plank_data, plank)
    end
    planks = NamedTuple{Tuple(plank_names)}(plank_data)
    return individuals(planks)
end

include("Advection/Advection.jl")
include("QuotaMode/QuotaMode.jl")
include("CarbonMode/CarbonMode.jl")
include("MacroMolecularMode/MacroMolecularMode.jl")
include("utils.jl")

using .Advection
import .Quota
import .Carbon
import .MacroMolecular

#####
##### some workarounds for function names
#####
construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN, mode::MacroMolecularMode) = 
    MacroMolecular.construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN)

construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN, mode::QuotaMode) = 
    Quota.construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN)

construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN, mode::CarbonMode) = 
    Carbon.construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN)

generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::MacroMolecularMode) =
    MacroMolecular.generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::QuotaMode) =
    Quota.generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::CarbonMode) =
    Carbon.generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture, mode::MacroMolecularMode) =
    MacroMolecular.plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture)
plankton_update!(plank, nuts, proc, p, plk, nothing::Nothing, ΔT, t, arch::Architecture, mode::MacroMolecularMode) =
    MacroMolecular.plankton_update!(plank, nuts, proc, p, plk, nothing::Nothing, ΔT, t, arch::Architecture)

plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture, mode::QuotaMode) =
    Quota.plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture)
plankton_update!(plank, nuts, proc, p, plk, nothing::Nothing, ΔT, t, arch::Architecture, mode::QuotaMode) =
    Quota.plankton_update!(plank, nuts, proc, p, plk, nothing::Nothing, ΔT, t, arch::Architecture)

plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture, mode::CarbonMode) =
    Carbon.plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture)
plankton_update!(plank, nuts, proc, p, plk, nothing::Nothing, ΔT, t, arch::Architecture, mode::CarbonMode) =
    Carbon.plankton_update!(plank, nuts, proc, p, plk, nothing::Nothing, ΔT, t, arch::Architecture)

end