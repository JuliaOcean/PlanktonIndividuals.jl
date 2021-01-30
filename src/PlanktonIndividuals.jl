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
include("$src"*"nut/fields.jl")
include("$src"*"nut/diffusivity.jl")
include("$src"*"nut/forcing.jl")
include("$src"*"nut/third_order_DSTFL.jl")
include("$src"*"nut/multi_dim_adv.jl")
include("$src"*"nut/halo_regions.jl")
include("$src"*"nut/gen_nut_fields.jl")
include("$src"*"nut/nutrient_processes.jl")
include("$src"*"plankton/gen_plankton.jl")
include("$src"*"plankton/advection/interpolation.jl")
include("$src"*"plankton/advection/advection_operations.jl")
include("$src"*"plankton/advection/plankton_advection.jl")
include("$src"*"plankton/advection/plankton_diffusion.jl")
include("$src"*"plankton/physiology/counts.jl")
include("$src"*"plankton/physiology/uptake_processes.jl")
include("$src"*"plankton/physiology/cell_division_loss.jl")
include("$src"*"plankton/physiology/plankton_operator.jl")
include("$src"*"plankton/physiology/plankton_update.jl")
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
