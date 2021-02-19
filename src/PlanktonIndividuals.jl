module PlanktonIndividuals

using NCDatasets, Serialization
using Random, Distributions, Statistics, Interpolations
using Printf, JLD2
using CUDA, KernelAbstractions, LinearAlgebra
using StructArrays

src=""

include("$src"*"architecture.jl")
include("$src"*"model/grids.jl")
include("$src"*"params/param_default.jl")
include("$src"*"params/param_update.jl")
include("$src"*"fields/fields.jl")
include("$src"*"fields/halo_regions.jl")
include("$src"*"diffusivity/tracer_diffusion.jl")
include("$src"*"diffusivity/plankton_diffusion.jl")
include("$src"*"advection/DST3FL.jl")
include("$src"*"advection/multi_dim_adv.jl")
include("$src"*"advection/interpolation.jl")
include("$src"*"advection/plankton_advection_kernels.jl")
include("$src"*"biogeochemistry/nutrient_forcings.jl")
include("$src"*"biogeochemistry/nutrient_fields.jl")
include("$src"*"biogeochemistry/nutrient_update.jl")
include("$src"*"plankton/plankton_generation.jl")
include("$src"*"plankton/plankton_advection.jl")
include("$src"*"plankton/physiology/counts.jl")
include("$src"*"plankton/physiology/uptake_processes.jl")
include("$src"*"plankton/physiology/division_loss_kernels.jl")
include("$src"*"plankton/physiology/division_loss.jl")
include("$src"*"plankton/physiology/plankton_growth.jl")
include("$src"*"plankton/plankton_physiology.jl")
include("$src"*"utils.jl")
include("$src"*"output/output_writers.jl")
include("$src"*"output/diagnostics.jl")
include("$src"*"model/timestepper.jl")
include("$src"*"model/models.jl")
include("$src"*"model/time_step.jl")
include("$src"*"model/simulations.jl")


export
    # model structures
    PI_Model, Grids, nutrient_fields,
    gen_Grid,
    Architecture, GPUs, CPUs,

    # read input functions
    read_IR_input, read_temp_input,
    update_bgc_params, update_phyt_params, 
    bgc_params_default, phyt_params_default,
    PrepRunDir,

    # initialize nutrient field and individual sets
    generate_nutrients,

    # Run the model
    PI_simulation, update!, vel_copy!, PI_TimeStep!,

    # write output functions
    write_nut_nc_each_step,
    write_individuals_to_jld2, write_diags_to_jld2,
    interior
end # module
