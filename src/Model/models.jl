mutable struct PlanktonModel
    arch::Architecture          # architecture on which models will run
    t::Float64                  # time in second
    iteration::Int64            # model interation
    individuals::individuals    # initial individuals generated by `setup_agents`
    nutrients::NamedTuple       # initial nutrient fields
    grid::AbstractGrid          # grid information
    bgc_params::Dict            # biogeochemical parameter set
    timestepper::timestepper    # operating Tuples and arrays for timestep
    mode::AbstractMode          # Carbon, Quota, or MacroMolecular
end

"""
    PlanktonModel(arch::Architecture, grid::AbstractGrid;
                  mode = QuotaMode(),
                  N_species = 1,
                  N_individual = [1024],
                  max_individuals = 1024*8,
                  bgc_params = nothing, 
                  phyt_params = nothing,
                  nut_initial = default_nut_init(),
                  t = 0.0,
                  )

Generate a `PlanktonModel` data structure. 

Keyword Arguments (Required)
============================
- `arch` : `CPU()` or `GPU()`. Computer architecture being used to run the model.
- `grid` : a `AbstractGrid` structure. Discrete grid for the model (resolution and geometry).

Keyword Arguments (Optional)
============================
- `mode` : Phytoplankton physiology mode, choose among CarbonMode(), QuotaMode(), or MacroMolecularMode().
- `N_species` : Number of species.
- `N_individual` : Number of individuals per species, should be a vector with `N_species` elements.
- `max_individuals` : Maximum number of individuals for each species the model can hold,
                    usually take the maximum of all the species and apply a factor to account for the growth
                    of individuals during one simulation.
- `bgc_params` : Parameter set for biogeochemical processes modeled in the model, use default if `nothing`, 
                    use `Dict` to update parameters, the format and names of parameters can be found by running `bgc_params_default()`.
- `phyt_params` : Parameter set for physiological processes of individuals modeled in the model, use default if `nothing`,
                    use `Dict` to update parameters, the format and names of parameters can be found by running `phyt_params_default(N_species, mode)`.
- `nut_initial` : The source of initial conditions of nutrient fields, should be either a `NamedTuple` 
                    or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
- `t` : Model time, start from 0 by default, in second.
"""
function PlanktonModel(arch::Architecture, grid::AbstractGrid;
                       mode = QuotaMode(),
                       N_species::Int64 = 1,
                       N_individual::Vector{Int64} = [1024],
                       max_individuals::Int64 = 8*1024,
                       bgc_params = nothing, 
                       phyt_params = nothing,
                       nut_initial = default_nut_init(),
                       t::Float64 = 0.0,
                       )

    @assert maximum(N_individual) ≤ max_individuals

    if arch == GPU() && !has_cuda()
        throw(ArgumentError("Cannot create a GPU model. No CUDA-enabled GPU was detected!"))
    end

    grid_d = replace_grid_storage(arch, grid)

    if isa(phyt_params, Nothing)
        phyt_params_final = phyt_params_default(N_species, mode)
    else
        if isa(phyt_params, Dict)
            phyt_params_final = update_phyt_params(phyt_params; N = N_species, mode = mode)
        else
            throw(ArgumentError("Phytoplankton parameters must be either Nothing or Dict!")) 
        end
    end

    if isa(bgc_params, Nothing)
        bgc_params_final = bgc_params_default()
    else
        if isa(bgc_params, Dict)
            bgc_params_final = update_bgc_params(bgc_params)
        else
            throw(ArgumentError("Phytoplankton parameters must be either Nothing or Dict!")) 
        end
    end

    inds = generate_individuals(phyt_params_final, arch, N_species, N_individual, max_individuals, grid_d, mode)

    nutrients = generate_nutrients(arch, grid_d, nut_initial)

    ts = timestepper(arch, grid_d, max_individuals)

    iteration  = 0

    model = PlanktonModel(arch, t, iteration, inds, nutrients, grid_d, bgc_params_final, ts, mode)

    return model
end

function show(io::IO, model::PlanktonModel)
    Nsp = length(model.individuals.phytos)
    N = Int(dot(model.individuals.phytos.sp1.data.ac,model.individuals.phytos.sp1.data.ac))
    cap = length(model.individuals.phytos.sp1.data.ac)
    print(io, "PlanktonModel:\n",
              "├── grid: $(short_show(model.grid))\n",
              "├── $(model.mode) is selected for phytoplankton physiology\n",
              "├── individuals: $(Nsp) phytoplankton species with $(N) individuals for each species\n",
              "└── maximum number of individuals: $(cap) per species")
end
