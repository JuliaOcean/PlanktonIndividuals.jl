module Model

export PlanktonModel
export PlanktonSimulation
export update!, vel_copy!

using CUDA
using StructArrays
using LinearAlgebra

using PlanktonIndividuals.Architectures
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields
using PlanktonIndividuals.Biogeochemistry
using PlanktonIndividuals.Parameters
using PlanktonIndividuals.Plankton
using PlanktonIndividuals.Output

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

import Base: show


include("utils.jl")
include("timestepper.jl")
include("models.jl")
include("time_step.jl")
include("simulations.jl")

end