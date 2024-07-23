module Model

export PlanktonModel
export TimeStep!

using StructArrays
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields
using PlanktonIndividuals.Biogeochemistry
using PlanktonIndividuals.Parameters
using PlanktonIndividuals.Plankton
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

import Base: show


include("timestepper.jl")
include("models.jl")
include("time_step.jl")

end