
using DataFrames, NetCDF, Printf, CSV, Serialization
using Random
using Distributions

#src="AgentPhytModel_3D/"
src=""

include("$src"*"parameters.jl")
include("$src"*"model_setup.jl")
include("$src"*"model_struct.jl")
include("$src"*"phyt_process.jl")
include("$src"*"utils.jl")
include("$src"*"agent_div.jl")
include("$src"*"dst3fl.jl")
include("$src"*"nutrient_processes.jl")
include("$src"*"2nd_adv_diffu.jl")

