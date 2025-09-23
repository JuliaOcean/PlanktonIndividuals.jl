module Individuals

export particle_advection!
export particle_diffusion!
export plankton_update!
export abiotic_particle_update!
export generate_individuals, individuals
export find_inds!, find_NPT!, acc_counts!, acc_chl!, calc_par!

export particle_interaction!, particle_release!
export particles_from_bcs!

using StructArrays
using Random
using KernelAbstractions

using PlanktonIndividuals.Architectures: device, Architecture, rng_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode
using PlanktonIndividuals: individuals, phytoplankton, abiotic_particle, phyto_setup, abiotic_setup, Palat

include("Advection/Advection.jl")
include("Plankton/QuotaMode/QuotaMode.jl")
include("Plankton/CarbonMode/CarbonMode.jl")
include("Plankton/MacroMolecularMode/MacroMolecularMode.jl")
include("Plankton/IronEnergyMode/IronEnergyMode.jl")
include("Abiotic/Abiotic.jl")
include("utils.jl")

using .Advection
using .Abiotic
import .Quota
import .Carbon
import .MacroMolecular
import .IronEnergy

#####
##### generate individuals of multiple species
#####
function generate_individuals(phyto::phyto_setup, abiotic::abiotic_setup, maxN, arch::Architecture, FT::DataType, g::AbstractGrid, mode::AbstractMode)
    plank_names = Symbol[]
    plank_data=[]

    for i in 1:phyto.Nsp
        name = Symbol("sp"*string(i))
        plank = construct_plankton(arch, i, phyto.params, maxN, FT, mode::AbstractMode)
        initialize_plankton!(plank, phyto.N[i], g, arch, mode)
        push!(plank_names, name)
        push!(plank_data, plank)
    end
    planks = NamedTuple{Tuple(plank_names)}(plank_data)

    ## add abiotic particles
    if isa(abiotic, Nothing)
        return individuals(planks, NamedTuple(;))
    else
        abiotic_names = Symbol[]
        abiotic_data = []

        for j in 1:abiotic.Nsa
            name = Symbol("sa"*string(j))
            particle =  construct_abiotic_particle(arch, j, abiotic.params, maxN, FT)
            initialize_abiotic_particle!(particle, abiotic.N[j], g, arch)
            push!(abiotic_names, name)
            push!(abiotic_data, particle)
        end
        abiotics = NamedTuple{Tuple(abiotic_names)}(abiotic_data)
        return individuals(planks, abiotics)
    end
end

#####
##### some workarounds for function names
#####
construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType, mode::MacroMolecularMode) = 
    MacroMolecular.construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)

construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType, mode::QuotaMode) = 
    Quota.construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)

construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType, mode::CarbonMode) = 
    Carbon.construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)

construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType, mode::IronEnergyMode) = 
    IronEnergy.construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)

initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::MacroMolecularMode) =
    MacroMolecular.initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::QuotaMode) =
    Quota.initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::CarbonMode) =
    Carbon.initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture, mode::IronEnergyMode) =
    IronEnergy.initialize_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)

plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::MacroMolecularMode) =
    MacroMolecular.plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::AbstractMode)

plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::QuotaMode) =
    Quota.plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::AbstractMode)

plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::CarbonMode) =
    Carbon.plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::AbstractMode)

plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::IronEnergyMode) =
    IronEnergy.plankton_update!(phyto, trs, proc, plk, diags_spcs, ΔT, t, arch::Architecture, mode::AbstractMode)

end
