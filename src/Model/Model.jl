module Model

import PlanktonIndividuals: CarbonMode, QuotaMode, MacroMolecularMode, Phytoplankton, phyto_setup,
    colony_setup, Palat, abiotic_setup, ModelOpts, timestepper, PlanktonModel, PlanktonDiagnostics

export PlanktonModel
export model_options

using StructArrays
using LinearAlgebra: dot

using PlanktonKernels.Architectures: CPU, GPU, Architecture, array_type, isfunctional
using PlanktonKernels.Grids: AbstractGrid, short_show, replace_grid_storage
using PlanktonKernels.Fields: Field, tracers_init, zero_fields!
using PlanktonKernels.Biogeochemistry: bgc_tracer_names, bgc_tracer_init, generate_bgc_tracers, bgc_tracer_update!, bgc_params_default, update_bgc_params

using PlanktonIndividuals.Individuals
using PlanktonIndividuals.Individuals: plankton_update!, colony_update!, generate_individuals,
    find_NPT!, acc_counts!, acc_chl!, calc_par!
using PlanktonIndividuals.Individuals.ParticleMotion: particle_motion!, colony_motion!
using PlanktonIndividuals.Individuals.Abiotic: particle_interaction!, particle_release!, particles_from_bcs!
using PlanktonIndividuals.Diagnostics
using PlanktonIndividuals.Diagnostics: diags_proc!


import Base: show


include("timestepper.jl")
include("models.jl")
include("time_step.jl")

end
