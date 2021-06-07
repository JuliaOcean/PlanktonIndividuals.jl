module PlanktonIndividuals

using NCDatasets, Serialization
using Random, Statistics, Interpolations
using Printf, JLD2
using CUDA, KernelAbstractions, CUDAKernels, Adapt
using LinearAlgebra
using StructArrays
using Pkg.Artifacts

p=dirname(pathof(PlanktonIndividuals))
artifact_toml = joinpath(p, "../Artifacts.toml")
surface_mixing_vels_hash = artifact_hash("surface_mixing_vels", artifact_toml)
surface_mixing_vels = joinpath(artifact_path(surface_mixing_vels_hash)*"/velocities.jld2")
default_temperature = joinpath(artifact_path(surface_mixing_vels_hash)*"/temp.bin")
default_PAR = joinpath(artifact_path(surface_mixing_vels_hash)*"/PAR.bin")

include("architecture.jl")
include("grids/regular_rectilinear_grid.jl")
include("grids/regular_lat_lon_grid.jl")
include("grids/vertically_stretched_lat_lon_grid.jl")
include("grids/areas_volumes_spacings.jl")
include("grids/utils.jl")
include("params/param_default.jl")
include("params/param_update.jl")
include("fields/boundary_conditions.jl")
include("fields/apply_bcs.jl")
include("fields/fields.jl")
include("fields/halo_regions.jl")
include("diffusivity/tracer_diffusion.jl")
include("diffusivity/plankton_diffusion.jl")
include("advection/DST3FL.jl")
include("advection/multi_dim_adv.jl")
include("advection/interpolation.jl")
include("plankton/plankton_advection_kernels.jl")
include("biogeochemistry/nutrient_forcings.jl")
include("biogeochemistry/nutrient_fields.jl")
include("biogeochemistry/nutrient_update.jl")
include("plankton/plankton_generation.jl")
include("plankton/plankton_advection.jl")
include("plankton/physiology/counts.jl")
include("plankton/physiology/uptake_processes.jl")
include("plankton/physiology/division_loss_kernels.jl")
include("plankton/physiology/division_loss.jl")
include("plankton/physiology/plankton_growth.jl")
include("plankton/plankton_physiology.jl")
include("utils.jl")
include("output/output_writers.jl")
include("output/diagnostics.jl")
include("output/diagnostics_plankton.jl")
include("model/timestepper.jl")
include("model/models.jl")
include("model/time_step.jl")
include("model/simulations.jl")


export
    # model structures
    PlanktonModel, RegularRectilinearGrid, RegularLatLonGrid,
    VerticallyStretchedLatLonGrid, LoadVerticallyStretchedLatLonGrid,
    Architecture, GPU, CPU, Periodic, Bounded,

    # read input functions
    read_IR_input, read_temp_input,
    update_bgc_params, update_phyt_params, 
    bgc_params_default, phyt_params_default,
    generate_nutrients,
    default_nut_init, PrepRunDir,

    # diagnostics
    PlanktonDiagnostics,

    # Run the model
    PlanktonSimulation, update!, vel_copy!, TimeStep!,

    # write output functions
    write_nut_nc_each_step,
    write_individuals_to_jld2, write_diags_to_jld2,
    interior
end # module
