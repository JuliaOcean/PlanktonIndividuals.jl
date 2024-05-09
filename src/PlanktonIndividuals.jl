module PlanktonIndividuals

if VERSION < v"1.6"
    error("This version of PlanktonIndividuals.jl requires Julia v1.6 or newer.")
end

export
    # Architectures
    Architecture, GPU, CPU,

    # Grids
    RectilinearGrid,
    LatLonGrid, LoadLatLonGrid,
    Periodic, Bounded,

    # Model
    PlanktonModel, 
    CarbonMode, QuotaMode, MacroMolecularMode,

    # BoundaryConditions
    set_bc!,

    # Simulation
    PlanktonSimulation, update!, vel_copy!, set_vels_fields!, set_PARF_fields!, set_temp_fields!,

    # Parameters
    default_PARF, default_temperature,
    update_bgc_params, update_phyt_params, 
    bgc_params_default, phyt_params_default,

    # Biogeochemistry
    generate_nutrients, default_nut_init,

    # Output
    PlanktonDiagnostics, PlanktonOutputWriter, interior,

    # Units
    second, minute, hour, meter, kilometer,
    seconds, minutes, hours, meters, kilometers,
    KiB, MiB, GiB, TiB

using CUDA
using Pkg.Artifacts

import Base: show

p=dirname(pathof(PlanktonIndividuals))
artifact_toml = joinpath(p, "../Artifacts.toml")
surface_mixing_vels_hash = artifact_hash("surface_mixing_vels", artifact_toml)
surface_mixing_vels = joinpath(artifact_path(surface_mixing_vels_hash)*"/velocities.jld2")
global_vels_hash = artifact_hash("OCCA_FlowFields", artifact_toml)
global_vels = joinpath(artifact_path(global_vels_hash)*"/OCCA_FlowFields.jld2")


"""
    AbstractMode
Abstract type for phytoplankton physiology modes supported by PlanktonIndividuals.
"""
abstract type AbstractMode end

"""
    CarbonMode <: AbstractMode
Type for the phytoplankton physiology mode which only resolves carbon quota.
"""
struct CarbonMode <: AbstractMode end

"""
    QuotaMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves carbon, nitrogen, and phosphorus quotas.
"""
struct QuotaMode <: AbstractMode end

"""
    MacroMolecularMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves marco-molecules.
"""
struct MacroMolecularMode <: AbstractMode end

include("Architectures.jl")
include("Units.jl")
include("Grids/Grids.jl")
include("Parameters/Parameters.jl")
include("Fields/Fields.jl")
include("Biogeochemistry/Biogeochemistry.jl")
include("Diagnostics/Diagnostics.jl")
include("Plankton/Plankton.jl")
include("Model/Model.jl")
include("Output/Output.jl")
include("Simulation/Simulation.jl")


using .Architectures
using .Grids
using .Parameters
using .Fields
using .Biogeochemistry
using .Diagnostics
using .Plankton
using .Model
using .Output
using .Simulation
using .Units

function __init__()
    if CUDA.has_cuda()
        @debug "CUDA-enabled GPU detected:"
        for (gpu, dev) in enumerate(CUDA.devices())
            @debug "$dev: $(CUDA.name(dev))"
        end
        CUDA.allowscalar(false)
    end
end

end # module
