module Model

export PlanktonModel
export TimeStep!

using CUDA
using StructArrays
using LinearAlgebra

using PlanktonIndividuals.Architectures
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields
using PlanktonIndividuals.Biogeochemistry
using PlanktonIndividuals.Parameters
using PlanktonIndividuals.Plankton
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

import Base: show


include("timestepper.jl")
include("models.jl")
include("time_step.jl")

end