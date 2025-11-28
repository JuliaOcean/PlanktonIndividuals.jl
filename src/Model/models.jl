mutable struct PlanktonModel
    arch::Architecture          # architecture on which models will run
    max_candidates::Int         # maximum number of candidate phytoplankton for interaction with one abiotic particle
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
                  max_individuals::Int = 8*1024,
                  max_abiotics::Int = 8*1024,
                  bgc_params = nothing, 
                  tracer_initial = default_tracer_init(),
                  phyto = nothing,
                  abiotic = nothing,
                  t::AbstractFloat = 0.0f0,
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
- `max_individuals` : Maximum number of individuals for each species the model can hold,
                    usually take the maximum of all the species and apply a factor to account for the growth
                    of individuals during one simulation.
- `max_abiotics` : Maximum number of abiotic particles for each species the model can hold.
- `bgc_params` : Parameter set for biogeochemical processes modeled in the model, use default if `nothing`, 
                    use `Dict` to update parameters, the format and names of parameters can be found by running `bgc_params_default()`.
- `tracer_initial` : The source of initial conditions of tracer fields, should be either a `NamedTuple` 
                    or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
- `phyto` : nothing or a `phyto_setup`. Whether to use default setup of phytoplankton in the model. If yes,
                    it should be a NamedTuple like this `phyto = phyto_setup(params = nothing, N = [2^10, 2^10], Nsp = 2)`.
- `abiotic` : nothing or a `abiotic_setup`. Whether to include abiotic particles in the model. If yes,
                    it should be a NamedTuple like this `abiotic = abiotic_setup(params = nothing, N = [2^10, 2^10], Nsa = 2, palat = [(:sp1, :sa1)])`.
- `t` : Model time, start from 0 by default, in second.
"""
function PlanktonModel(arch::Architecture, grid::AbstractGrid;
                       FT = Float32,
                       mode = QuotaMode(),
                       max_individuals::Int = 8*1024,
                       max_abiotics::Int = 8*1024,
                       bgc_params = nothing, 
                       tracer_initial = default_tracer_init(),
                       phyto = nothing,
                       abiotic = nothing,
                       t::AbstractFloat = 0.0f0,
                       max_candidates::Int = 25,
                       )

    @assert isfunctional(arch) == true

    if isa(bgc_params, Nothing)
        bgc_params = bgc_params_default(FT)
        bgc_params_final = update_bgc_params(bgc_params, FT)
    elseif isa(bgc_params, Dict)
        bgc_params_final = update_bgc_params(bgc_params, FT)
    else
        throw(ArgumentError("Phytoplankton parameters must be either Nothing or Dict!")) 
    end

    grid_d = replace_grid_storage(arch, grid)

    if isa(phyto, Nothing)
        N_inds = [2^10]
        Nsp = 1
        phyt_params = phyt_params_default(Nsp, mode)
        phyt_params = update_phyt_params(phyt_params, FT; N = Nsp, mode = mode)
        phyto = phyto_setup(phyt_params, N_inds, Nsp)
    elseif isa(phyto, phyto_setup)
        @assert maximum(phyto.N) ≤ max_individuals
        if length(phyto.N) ≠ phyto.Nsp
            throw(ArgumentError("PlanktonModel: `phyto`: The length of `N` must be $(phyto.Nsp), the same as `Nsp`, each species has its own initial condition"))
        end
        if isa(phyto.params, Nothing)
            phyto.params = phyt_params_default(phyto.Nsp, mode)
            phyto.params = update_phyt_params(phyto.params, FT; N = phyto.Nsp, mode = mode)
        elseif isa(phyto.params, Dict)
            phyto.params = update_phyt_params(phyto.params, FT; N = phyto.Nsp, mode = mode)
        else
            throw(ArgumentError("Phytoplankton parameters must be either Nothing or Dict!")) 
        end
    else
        throw(ArgumentError("PlanktonModel:`phyto` must be either Nothing or `phyto_setup`!")) 
    end

    if isa(abiotic, Nothing)
        intac = nothing
    elseif isa(abiotic, abiotic_setup)
        @assert maximum(abiotic.N) ≤ max_abiiotics
        intac = zeros(Int, max_candidates, max_abiotics) |> array_type(arch)
        if length(abiotic.N) ≠ abiotic.Nsa
            throw(ArgumentError("PlanktonModel: `abiotic`: The length of `N` must be $(abiotic.Nsa), the same as `Nsa`, each species has its own initial condition"))
        end
        if isa(abiotic.params, Nothing)
            abiotic.params = abiotic_params_default(abiotic.Nsa)
            abiotic.params = update_abiotic_params(abiotic.params, FT; N = abiotic.Nsa)
        elseif isa(abiotic.params, Dict)
            abiotic.params = update_abiotic_params(abiotic.params, FT; N = abiotic.Nsa)
        else
            throw(ArgumentError("Abiotic particle parameters must be either Nothing or Dict!")) 
        end
    else
        throw(ArgumentError("PlanktonModel:`abiotic` must be either Nothing or `abiotic_setup`!")) 
    end

    inds = generate_individuals(phyto, abiotic, max_individuals, max_abiotics, arch, FT, grid_d, mode)

    ##### check palatability between plank and abiotic
    if isa(abiotic, Nothing)
        palat = Palat([], [])
    else
        SPs = keys(inds.phytos)
        SAs = keys(inds.abiotics)
        for p in abiotic.palat.intac
            if p[1] ∉ SPs
                throw(ArgumentError("Abiotic: $(p[1]) is not generated"))
            end
            if p[2] ∉ SAs
                throw(ArgumentError("Abiotic: $(p[s]) is not generated"))
            end
        end
        for p in abiotic.palat.release
            if p[1] ∉ SPs
                throw(ArgumentError("Abiotic: $(p[1]) is not generated"))
            end
            if p[2] ∉ SAs
                throw(ArgumentError("Abiotic: $(p[s]) is not generated"))
            end
        end
        palat = abiotic.palat
    end

    tracers = generate_tracers(arch, grid_d, tracer_initial, FT)

    ts = timestepper(arch, FT, grid_d, max_individuals,max_abiotics, intac, palat)

    iteration  = 0

    model = PlanktonModel(arch, max_candidates, FT, t, iteration, inds, tracers, grid_d, bgc_params_final, ts, mode)

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