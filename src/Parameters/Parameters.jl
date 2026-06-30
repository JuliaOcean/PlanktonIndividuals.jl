module Parameters

export bgc_params_default
export phyt_params_default, colony_params_default
export abiotic_params_default
export default_PARF, default_temperature
export update_bgc_params
export update_phyt_params, update_colony_params
export update_abiotic_params

using PlanktonIndividuals.Grids

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

include("param_default.jl")
include("param_update.jl")

end