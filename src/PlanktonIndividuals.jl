module PlanktonIndividuals

if VERSION < v"1.11"
    error("This version of PlanktonIndividuals.jl requires Julia v1.11 or newer.")
end

export
    # Types and Structs
    CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode, ProteinMode,
    phyto_setup, colony_setup, abiotic_setup, Palat,

    # Model
    PlanktonModel, model_options,

    # BoundaryConditions
    set_bc_particle!,

    # Simulation
    PlanktonSimulation, update!, set_vels_fields!, set_PARF_fields!, set_temp_fields!,

    # Mode defaults and model parameter updates
    default_PARF, default_temperature,
    update_phyt_params, phyt_params_default,
    update_abiotic_params, colony_params_default,
    update_colony_params, abiotic_params_default,

    # Output
    PlanktonDiagnostics, PlanktonOutputWriter

using Artifacts: artifact_hash, artifact_path


p=dirname(pathof(PlanktonIndividuals))
artifact_toml = joinpath(p, "../Artifacts.toml")
surface_mixing_vels_hash = artifact_hash("surface_mixing_vels", artifact_toml)
surface_mixing_vels = joinpath(artifact_path(surface_mixing_vels_hash)*"/velocities.jld2")
global_vels_hash = artifact_hash("OCCA_FlowFields", artifact_toml)
global_vels = joinpath(artifact_path(global_vels_hash)*"/OCCA_FlowFields.jld2")


using PlanktonKernels.Architectures: Architecture
using PlanktonKernels.Grids: AbstractGrid
using PlanktonKernels.Fields: BoundaryConditions

include("model_structs.jl")
include("Diagnostics/Diagnostics.jl")
include("Individuals/Individuals.jl")
include("Model/Model.jl")
include("Output/Output.jl")
include("Simulation/Simulation.jl")

using .Individuals.Abiotic: set_bc_particle!

using .Diagnostics
using .Individuals
using .Individuals: IndividualParticles
using .Model
using .Output
using .Simulation

end # module
