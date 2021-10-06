module PlanktonIndividuals

if VERSION < v"1.6"
    error("This version of PlanktonIndividuals.jl requires Julia v1.6 or newer.")
end

export
    # Model
    PlanktonModel, PlanktonSimulation, update!, vel_copy!,
    CarbonMode, QuotaMode, MacroMolecularMode,

    # Grids
    RegularRectilinearGrid, RegularLatLonGrid,
    VerticallyStretchedLatLonGrid, LoadVerticallyStretchedLatLonGrid,
    Periodic, Bounded,

    # Architectures
    Architecture, GPU, CPU,

    # Parameters
    default_PARF, default_temperature,
    update_bgc_params, update_phyt_params, 
    bgc_params_default, phyt_params_default,

    # Biogeochemistry
    generate_nutrients, default_nut_init,

    # Output
    PlanktonDiagnostics, PrepRunDir, interior


using CUDA
using Pkg.Artifacts

import Base: show

p=dirname(pathof(PlanktonIndividuals))
artifact_toml = joinpath(p, "../Artifacts.toml")
surface_mixing_vels_hash = artifact_hash("surface_mixing_vels", artifact_toml)
surface_mixing_vels = joinpath(artifact_path(surface_mixing_vels_hash)*"/velocities.jld2")

abstract type AbstractMode end
struct CarbonMode <: AbstractMode end
struct QuotaMode <: AbstractMode end
struct MacroMolecularMode <: AbstractMode end

include("Architectures.jl")
include("Grids/Grids.jl")
include("Parameters/Parameters.jl")
include("Fields/Fields.jl")
include("Biogeochemistry/Biogeochemistry.jl")
include("Output/Output.jl")
include("Plankton/Plankton.jl")
include("Model/Model.jl")


using .Architectures
using .Grids
using .Parameters
using .Fields
using .Biogeochemistry
using .Output
using .Plankton
using .Model

function __init__()
    if CUDA.has_cuda()
        CUDA.allowscalar(false)
    end
end

end # module