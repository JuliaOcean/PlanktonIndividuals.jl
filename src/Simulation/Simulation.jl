module Simulation

export PlanktonSimulation
export update!, vel_copy!, set_vels_fields!, set_PARF_fields!, set_temp_fields!

using StructArrays
using LinearAlgebra: dot

using PlanktonIndividuals.Architectures
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Fields
using PlanktonIndividuals.Biogeochemistry
using PlanktonIndividuals.Parameters
using PlanktonIndividuals.Plankton
using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Model
using PlanktonIndividuals.Output

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode, IronEnergyMode

import Base: show


include("simulations.jl")
include("utils.jl")
include("update.jl")

end