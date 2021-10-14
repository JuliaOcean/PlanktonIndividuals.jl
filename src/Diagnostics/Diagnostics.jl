module Diagnostics

export PlanktonDiagnostics
export diags_spcs!, diags_proc!
export tracer_avail_diags, plank_avail_diags

using CUDA

using PlanktonIndividuals.Architectures
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

import Base: show


include("diagnostics_struct.jl")
include("diagnostics_kernels.jl")

end