module PlanktonIndividuals

using NetCDF, Serialization
using Random, Distributions, Interpolations, Statistics
using Printf

src=""

include("$src"*"model_struct.jl")
include("$src"*"model_setup.jl")
include("$src"*"phyt_process.jl")
include("$src"*"zooplankton.jl")
include("$src"*"utils.jl")
include("$src"*"output_writers.jl")
include("$src"*"agent_div.jl")
include("$src"*"dst3fl.jl")
include("$src"*"nutrient_processes.jl")
include("$src"*"2nd_adv_diffu.jl")
include("$src"*"models.jl")
include("$src"*"time_step.jl")
include("$src"*"param_default.jl")


export
    # model structures
    PI_Model, grids, nutrient_fields, velocity,
    RunOptions, RunParams, read_Ogrids,

    # read input functions
    read_IR_input, read_temp_input,
    update_params!, grid_offline, param_default,
    PrepRunDir, generate_vel_itp,

    # initialize nutrient field and individual sets
    setup_agents, setup_nutrients, load_nut_initials,

    # Run the model
    RunParam, RunOption, PI_TimeStep!,
    PI_advectRK4!, PI_advect!,

    # write output functions
    write_nut_nc_alltime, write_nut_nc_each_step,
    count_vertical_num, count_horizontal_num,
    write_output
end # module
