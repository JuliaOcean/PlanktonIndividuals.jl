module Output

export PlanktonDiagnostics
export diags_spcs!, diags_proc!
export PrepRunDir
export write_diags_to_jld2, write_individuals_to_jld2
export write_species_dynamics
export tracer_avail_diags, plank_avail_diags

using CUDA
using LinearAlgebra, Statistics, JLD2, NCDatasets, Serialization, Printf

using PlanktonIndividuals.Architectures
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

import Base: show


include("output_writers.jl")
include("diagnostics.jl")
include("diagnostics_plankton.jl")

end