mutable struct PlanktonModel
    arch::Architecture          # architecture on which models will run
    FT::DataType                # floating point data type
    t::AbstractFloat            # time in second
    iteration::Int              # model interation
    individuals::individuals    # individuals
    tracers::NamedTuple         # tracer fields
    grid::AbstractGrid          # grid information
    bgc_params::Dict            # biogeochemical parameter set
    timestepper::timestepper    # operating Tuples and arrays for timestep
    mode::AbstractMode          # Carbon, Quota, or MacroMolecular
end

"""
    PlanktonModel(arch::Architecture, grid::AbstractGrid;
                  FT = Float32,
                  mode = QuotaMode(),
                  N_species = 1,
                  N_individual = [1024],
                  max_individuals = 1024*8,
                  bgc_params = nothing, 
                  phyt_params = nothing,
                  tracer_initial = default_tracer_init(),
                  abiotic = false,
                  t = 0.0,
                  )

Generate a `PlanktonModel` data structure. 

Keyword Arguments (Required)
============================
- `arch` : `CPU()` or `GPU()`. Computer architecture being used to run the model.
- `grid` : a `AbstractGrid` structure. Discrete grid for the model (resolution and geometry).

Keyword Arguments (Optional)
============================
- `FT`: Floating point data type. Default: `Float32`.
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
- `tracer_initial` : The source of initial conditions of tracer fields, should be either a `NamedTuple` 
                    or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
- `abiotic` : false or a NamedTuple. Whether to include abiotic particles in the model. If yes,
                    it should be a NamedTuple like this `abiotic = (params = nothing, N = [2^10, 2^10], Nsa = 2, palat = [(:sp1, :sa1)])`.
- `t` : Model time, start from 0 by default, in second.
"""
function PlanktonModel(arch::Architecture, grid::AbstractGrid;
                       FT = Float32,
                       mode = QuotaMode(),
                       N_species::Int = 1,
                       N_individual::Vector{Int} = [1024],
                       max_individuals::Int = 8*1024,
                       bgc_params = nothing, 
                       phyt_params = nothing,
                       tracer_initial = default_tracer_init(),
                       abiotic = nothing,
                       t::AbstractFloat = 0.0f0,
                       )

    @assert maximum(N_individual) ≤ max_individuals

    @assert isfunctional(arch) == true

    grid_d = replace_grid_storage(arch, grid)

    if isa(phyt_params, Nothing)
        phyt_params = phyt_params_default(N_species, mode)
        phyt_params_final = update_phyt_params(phyt_params, FT; N = N_species, mode = mode)
    elseif isa(phyt_params, Dict)
        phyt_params_final = update_phyt_params(phyt_params, FT; N = N_species, mode = mode)
    else
        throw(ArgumentError("Phytoplankton parameters must be either Nothing or Dict!")) 
    end

    if isa(bgc_params, Nothing)
        bgc_params = bgc_params_default(FT)
        bgc_params_final = update_bgc_params(bgc_params, FT)
    elseif isa(bgc_params, Dict)
        bgc_params_final = update_bgc_params(bgc_params, FT)
    else
        throw(ArgumentError("Phytoplankton parameters must be either Nothing or Dict!")) 
    end

    if isa(abiotic, NamedTuple)
        if isa(abiotic.params, Nothing)
            abiotic_params = abiotic_params_default(abiotic.Nsa)
            abiotic_params_final = update_abiotic_params(abiotic_params, FT; N = abiotic.Nsa)
        elseif isa(abiotic.params, Dict)
            abiotic_params_final = update_abiotic_params(abiotic_params, FT; N = abiotic.Nsa)
        else
            throw(ArgumentError("Abiotic particle parameters must be either Nothing or Dict!")) 
        end
        abiotic_final = (params = abiotic_params_final, Nsa = abiotic.Nsa, N = abiotic.N, palat = abiotic.palat)
    else
        abiotic_final = nothing
    end

    inds = generate_individuals(phyt_params_final, arch, N_species, N_individual, max_individuals, FT, grid_d, mode; abiotic = abiotic_final)

    ##### check palatability between plank and abiotic
    if isa(abiotic_final, Nothing)
        palat = []
    else
        SPs = keys(inds.phytos)
        SAs = keys(inds.abiotics)
        for p in abiotic_final.palat
            if p[1] ∉ SPs
                throw(ArgumentError("Abiotic: $(p[1]) is not generated"))
            end
            if p[2] ∉ SAs
                throw(ArgumentError("Abiotic: $(p[s]) is not generated"))
            end
        end
    end

    tracers = generate_tracers(arch, grid_d, tracer_initial, FT)

    ts = timestepper(arch, FT, grid_d, max_individuals, abiotic_final.palat)

    iteration  = 0

    model = PlanktonModel(arch, FT, t, iteration, inds, tracers, grid_d, bgc_params_final, ts, mode)

    return model
end

function calc_active_individuals(particle)
    Nsp = length(particle)
    N = zeros(Int, Nsp)
    for i in 1:Nsp
        N[i] = Int(dot(particle[i].data.ac, particle[i].data.ac))
    end
    return N
end

function show(io::IO, model::PlanktonModel)
    Nsp = length(model.individuals.phytos)
    N = calc_active_individuals(model.individuals.phytos)
    cap = length(model.individuals.phytos.sp1.data.ac)
    if model.individuals.abiotics == NamedTuple(;)
        s = "├── No abiotic particles available"
    else
        abiotic_Nsp = length(model.individuals.abiotics)
        abiotic_N = calc_active_individuals(model.individuals.abiotics)
        s = "├── a biotic particles: $(abiotic_Nsp) species with $(abiotic_N) individuals for each species\n"
    end

    print(io, "PlanktonModel:\n",
              "├── floating point data type: $(model.FT)\n",
              "├── grid: $(short_show(model.grid))\n",
              "├── $(model.mode) is selected for phytoplankton physiology\n",
              "├── phytoplankton: $(Nsp) species with $(N) individuals for each species\n",
              s,
              "└── maximum number of individuals: $(cap) per species")
end