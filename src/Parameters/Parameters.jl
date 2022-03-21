module Parameters

export bgc_params_default, phyt_params_default
export default_PARF, default_temperature
export update_bgc_params, update_phyt_params

using CUDA
using PlanktonIndividuals.Grids

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

include("param_default.jl")
include("param_update.jl")

end