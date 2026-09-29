module Diagnostics

export PlanktonDiagnostics
export diags_spcs!, diags_proc!, diags_colony!

using KernelAbstractions
using StructArrays

using PlanktonIndividuals.Architectures: device, Architecture, array_type
using PlanktonIndividuals.Grids

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode
using PlanktonIndividuals: individuals, phytoplankton, abiotic_particle

import Base: show

include("diagnostics_struct.jl")
include("diagnostics_kernels.jl")

end