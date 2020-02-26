module PhytoAgentModel

greet() = print("Hello World!")

using DataFrames, NetCDF, JLD, Serialization, CSV
using Random, Distributions, Interpolations
using Printf, YAML

src=""

include("$src"*"model_struct.jl")
include("$src"*"model_setup.jl")
include("$src"*"phyt_process.jl")
include("$src"*"zooplankton.jl")
include("$src"*"utils.jl")
include("$src"*"agent_div.jl")
include("$src"*"dst3fl.jl")
include("$src"*"nutrient_processes.jl")
include("$src"*"2nd_adv_diffu.jl")
include("$src"*"model.jl")
include("$src"*"option_params.jl")


export
    # model structures
    PA_Model, grids, nutrient_fields, velocity,
    RunOptions, RunParams, PlankOpt, read_Ogrids,

    # read input functions
    read_default_IR_input, read_default_temp_input,
    update_params, grid_offline, param_default,
    PrepRunDir, generate_vel_itp,

    # initialize nutrient field and individual sets
    setup_agents, setup_nutrients,

    # Run the model
    PA_ModelRun, RunParam, RunOption, PA_TimeStep!,
    PA_advectRK4!, PA_advect!,

    # write output functions
    create_output, sort_species, convert_coordinates,
    count_vertical_num, count_horizontal_num,
    compute_mean_species, testB1B2
end # module
