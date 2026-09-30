module Simulation

import PlanktonIndividuals: PlanktonInput, PlanktonSimulation

export PlanktonSimulation
export default_PARF, default_temperature
export update!, set_vels_fields!, set_PARF_fields!, set_temp_fields!

using PlanktonKernels.Fields: vel_copy!, copy_interior!, validate_bcs

using StructArrays
using LinearAlgebra: dot

using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids: AbstractGrid, Bounded, replace_grid_storage
import PlanktonKernels.Grids: short_show

using PlanktonIndividuals.Individuals
using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Model
using PlanktonIndividuals.Model: TimeStep!
using PlanktonIndividuals.Output
using PlanktonIndividuals.Output: write_output!, humanize_filesize


import Base: show


include("default_forcing.jl")
include("simulations.jl")
include("utils.jl")
include("update.jl")

end
