mutable struct PlanktonModel
    arch::Architecture          # architecture on which models will run
    t::Int64                    # time in second
    iteration::Int64            # model interation
    individuals::individuals    # initial individuals generated by `setup_agents`
    nutrients::NamedTuple       # initial nutrient fields
    grid::AbstractGrid          # grid information
    bgc_params::Dict            # biogeochemical parameter set
    timestepper::timestepper    # operating Tuples and arrays for timestep
end

"""
    PlanktonModel(arch::Architecture, grid::RegularRectilinearGrid;
                  N_species = 1,
                  N_individual = 1024,
                  max_individuals = 8*1024,
                  bgc_params = bgc_params_default(), 
                  phyt_params = phyt_params_default(),
                  nut_initial = default_nut_init(),
                  t = 0.0,
                  mask = nothing,
                  )

Generate a `PlanktonModel` data structure. 

Keyword Arguments (Required)
============================
- `arch` : `CPU()` or `GPU()`. Computer architecture being used to run the model.
- `grid` : a `AbstractGrid` structure. Discrete grid for the model (resolution and geometry).

Keyword Arguments (Optional)
============================
- `N_species` : Number of species.
- `N_individual` : Number of individuals per species.
- `max_individuals` : Maximum number of individuals per species the model can hold.
- `bgc_params` : Parameter set for biogeochemical processes modeled in the model.
- `phyt_params` : Parameter set for physiological processes of individuals modeled in the model.
- `nut_initial` : The source of initial conditions of nutrient fields, should be either a `NamedTuple` 
                           or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
- `t` : Model time, start from 0 by default, in second.
- `mask` : Mask out the individuals and tracers generated out of the domain, a 3D array with size `(Nx, Ny, Nz)`.
"""
function PlanktonModel(arch::Architecture, grid::RegularRectilinearGrid;
                       N_species::Int64 = 1,
                       N_individual::Int64 = 1024,
                       max_individuals::Int64 = 8*1024,
                       bgc_params = bgc_params_default(), 
                       phyt_params = phyt_params_default(),
                       nut_initial = default_nut_init(),
                       t::Int64 = 0,
                       mask = nothing,
                       )

    if arch == GPU() && !has_cuda()
        throw(ArgumentError("Cannot create a GPU model. No CUDA-enabled GPU was detected!"))
    end

    inds = individuals(phyt_params, arch, N_species, N_individual, max_individuals)

    for plank in inds.phytos
        gen_individuals!(plank, N_individual, grid, arch; mask = mask)
    end

    nutrients = generate_nutrients(arch, grid, nut_initial; mask = mask)

    ts = timestepper(arch, grid, N_individual, max_individuals)

    iteration  = 0

    model = PlanktonModel(arch, t, iteration, inds, nutrients, grid, bgc_params, ts)

    return model
end

import Base: show

function show(io::IO, model::PlanktonModel)
    Nsp = length(model.individuals.phytos)
    N = Int(dot(model.individuals.phytos.sp1.data.ac,model.individuals.phytos.sp1.data.ac))
    cap = length(model.individuals.phytos.sp1.data.ac)
    print(io, "grid: Nx = $(model.grid.Nx), Ny = $(model.grid.Ny), Nz = $(model.grid.Nz)\n",
              "individuals: $(Nsp) phytoplankton species each with $(N) individuals\n",
              "maximum number of individuals: $(cap) per species\n")
end