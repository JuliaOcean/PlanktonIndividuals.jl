module Simulation

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
using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Model
using PlanktonIndividuals.Output

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

import Base: show


include("utils.jl")
include("simulations.jl")
include("update.jl")

end