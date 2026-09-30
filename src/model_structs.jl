# Shared type definitions, ordered so field types are defined before use.

"""
    AbstractMode
Abstract type for phytoplankton physiology modes supported by PlanktonIndividuals.
"""
abstract type AbstractMode end

"""
    CarbonMode <: AbstractMode
Type for the phytoplankton physiology mode which only resolves carbon quota.
"""
struct CarbonMode <: AbstractMode end

"""
    QuotaMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves carbon, nitrogen, and phosphorus quotas.
"""
struct QuotaMode <: AbstractMode end

"""
    MacroMolecularMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves marco-molecules.
"""
struct MacroMolecularMode <: AbstractMode end

"""
    IronEnergyMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves carbon, nitrogen, phosphorus, and iron quotas. This mode also resolves energy.
"""
struct IronEnergyMode <: AbstractMode end

"""
    ProteinMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves protein synthesis.
"""
struct ProteinMode <: AbstractMode end

##### struct for phytoplankton
mutable struct Phytoplankton
    data::AbstractArray
    p::NamedTuple
end

##### struct for colony
mutable struct ColonyParticle
    spcs::NamedTuple
    intac::AbstractArray
end

##### struct for abiotic particles
mutable struct AbioticParticle
    data::AbstractArray
    p::NamedTuple
    bc::BoundaryConditions
end

struct IndividualParticles
    phytos::NamedTuple
    abiotics::NamedTuple
    colonies::NamedTuple
end

mutable struct phyto_setup
    params::Union{Nothing, Dict}
    N::AbstractArray
    Nsp::Int64
end

mutable struct colony_setup
    params::Union{Nothing, AbstractArray}
    N::AbstractArray
    Nsp::AbstractArray
    Ncl::Int64
end

mutable struct Palat
    intac::AbstractArray
    release::AbstractArray
end
mutable struct abiotic_setup
    params::Union{Nothing, Dict}
    N::AbstractArray
    Nsa::Int64
    palat::Palat
end

mutable struct timestepper
    Gcs::NamedTuple                         # a NamedTuple same as tracers to store tendencies
    tracer_temp::NamedTuple                 # a NamedTuple same as tracers to store tracers fields in multi-dims advection scheme
    vel₀::NamedTuple                        # a NamedTuple with u, v, w velocities
    vel½::NamedTuple                        # a NamedTuple with u, v, w velocities
    vel₁::NamedTuple                        # a NamedTuple with u, v, w velocities
    PARF::AbstractArray                     # a (Cu)Array to store surface PAR field of each timestep
    temp::AbstractArray                     # a (Cu)Array to store temperature field of each timestep
    flux_sink::AbstractArray                # a (Cu)Array to store sinking flux field of each timestep
    plk::NamedTuple                         # a NamedTuple same as tracers to store interactions with individuals
    par::AbstractArray                      # a (Cu)Array to store PAR field of each timestep
    par₀::AbstractArray                     # a (Cu)Array to store PAR field of the previous timestep
    Chl::AbstractArray                      # a (Cu)Array to store Chl field of each timestep
    pop::AbstractArray                      # a (Cu)Array to store population field of each timestep
    rnd::AbstractArray                      # a StructArray of random numbers for plankton diffusion or grazing, mortality and division.
    rnd_3d::AbstractArray                   # a (Cu)Array of random numbers for tracer-particle interaction
    velos::AbstractArray                    # a StructArray of intermediate values for RK4 particle advection
    trs::AbstractArray                      # a StructArray of tracers of each individual
    intac::Union{Nothing, AbstractArray}    # Top-K candidate phyto IDs for abiotic particles
    palat::Palat                            # a `Palat` to store the interaction between species
end

mutable struct ModelOpts
    max_individuals::Int        # maximum number of individuals for each species the model can hold
    max_candidates::Int         # maximum number of candidate phytoplankton for interaction with one abiotic particle
    kc::Float64                 # light attenuation coefficient, unit: m^-1
    kw::Float64                 # light attenuation coefficient, unit: m^-1
    shared_graz::Float64        # whether to use shared grazing for all species, 1: yes, 0: no
end

mutable struct PlanktonModel
    arch::Architecture                   # architecture on which models will run
    options::ModelOpts                   # model options
    FT::DataType                         # floating point data type
    t::AbstractFloat                     # time in second
    iteration::Int                       # model interation
    individuals::IndividualParticles     # individuals
    tracers::NamedTuple                  # tracer fields
    grid::AbstractGrid                   # grid information
    bgc_params::Dict                     # biogeochemical parameter set
    timestepper::timestepper             # operating Tuples and arrays for timestep
    mode::AbstractMode                   # Carbon, Quota, or MacroMolecular
end

mutable struct PlanktonDiagnostics
    phytos::NamedTuple       # for each species of phytoplankton
    abiotics::NamedTuple     # for each species of abiotic particle
    colonies::NamedTuple     # for each colony
    tracer::NamedTuple       # for tracers
    iteration_interval::Int  # time interval that the diagnostics is time averaged
end

mutable struct PlanktonOutputWriter
    filepath::String
    write_log::Bool
    save_diags::Bool
    save_phytoplankton::Bool
    save_abiotic_particle::Bool
    save_colony::Bool
    diags_file::String
    phytoplankton_file::String
    phytoplankton_include::Tuple
    phytoplankton_iteration_interval::Int
    abiotic_particle_file::String
    abiotic_particle_include::Tuple
    abiotic_particle_iteration_interval::Int
    colony_file::String
    colony_include::Tuple
    colony_iteration_interval::Int
    max_filesize::Number # in Bytes
    part_diags::Int
    part_phytoplankton::Int
    part_abiotic_particle::Int
    part_colony::Int
end

mutable struct PlanktonInput
    temp::AbstractArray{AbstractFloat,4}      # temperature
    PARF::AbstractArray{AbstractFloat,3}      # PARF
    vels::NamedTuple                          # velocity fields
    ΔT_vel::AbstractFloat                     # time step of velocities provided
    ΔT_PAR::AbstractFloat                     # time step of surface PAR provided
    ΔT_temp::AbstractFloat                    # time step of temperature provided
end

mutable struct PlanktonSimulation
    model::PlanktonModel                                # Model object
    input::PlanktonInput                                # model input, temp, PAR, and velocities
    diags::Union{PlanktonDiagnostics,Nothing}           # diagnostics
    ΔT::AbstractFloat                                   # model time step
    iterations::Int                                     # run the simulation for this number of iterations
    output_writer::Union{PlanktonOutputWriter,Nothing}  # Output writer
end
