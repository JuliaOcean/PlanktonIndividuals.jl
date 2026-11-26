module Output

import PlanktonIndividuals: PlanktonOutputWriter

export PlanktonOutputWriter

using LinearAlgebra: dot
using JLD2
using Printf: @sprintf

using PlanktonKernels.Fields: interior

using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Model
using PlanktonIndividuals: CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

import Base: show


include("output_writers.jl")
include("write_outputs.jl")

end
