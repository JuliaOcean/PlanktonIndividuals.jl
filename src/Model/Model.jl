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
using PlanktonIndividuals.Individuals
using PlanktonIndividuals.Diagnostics

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode
using PlanktonIndividuals: individuals, phytoplankton, abiotic_particle

import Base: show


include("timestepper.jl")
include("models.jl")
include("time_step.jl")

end