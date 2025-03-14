module Parameters

export bgc_params_default, phyt_params_default, abiotic_params_default
export default_PARF, default_temperature
export update_bgc_params, update_phyt_params, update_abiotic_params

using PlanktonIndividuals.Grids

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

include("param_default.jl")
include("param_update.jl")

end