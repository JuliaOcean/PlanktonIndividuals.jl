mutable struct PlanktonDiagnostics
    phytos::NamedTuple       # for each species of phytoplankton
    abiotics::NamedTuple     # for each species of abiotic particle
    colonies::NamedTuple     # for each colony
    tracer::NamedTuple       # for tracers
    iteration_interval::Int  # time interval that the diagnostics is time averaged
end

"""
    PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                               phytoplankton = (:num, :graz, :mort, :dvid, :ptc),
                               abiotic_particle = (),
                               colony = (),
                               time_interval = 1)

Generate a `PlanktonDiagnostics` structure.

Keyword Arguments (Optional)
============================
- `tracer` : a `Tuple` containing the names of nutrient fields to be diagnosed.
- `phytoplankton` : a `Tuple` containing the names of physiological processes of phytoplankton individuals to be diagnosed.
- `abiotic_particle` : a `Tuple` containing the names of state variables of abiotic particles to be diagnosed.
- `colony` : a `Tuple` containing the names of physiological processes of colony to be diagnosed.
- `iteration_interval` : The number of timesteps that diagnostics is averaged, 1 iteration by default.
"""
function PlanktonDiagnostics(model; tracer=(),
                                    phytoplankton = (:num, :graz, :mort, :dvid, :ptc),
                                    abiotic_particle = (),
                                    colony = (),
                                    iteration_interval::Int = 1)
    
    @assert isa(tracer, Tuple)
    @assert isa(phytoplankton, Tuple)
    @assert isa(abiotic_particle, Tuple)
    @assert isa(colony, Tuple)

    diag_avail(tracer, phytoplankton, abiotic_particle, colony, model)

    ntr   = length(tracer)
    nproc_phyto = length(phytoplankton)
    nproc_abiotic = length(abiotic_particle)
    nproc_colony = length(colony)

    trs   = []
    phyto_procs = []
    abiotic_procs = []
    colony_procs = []
    FT = model.FT
    total_size = (model.grid.Nx+model.grid.Hx*2, model.grid.Ny+model.grid.Hy*2, model.grid.Nz+model.grid.Hz*2)

    ##### tracers
    for i in 1:ntr
        tr = zeros(FT, total_size) |> array_type(model.arch)
        push!(trs, tr)
    end
    tr_d1 = zeros(FT, total_size) |> array_type(model.arch)
    tr_default = (PAR = tr_d1,)
    diag_tr = NamedTuple{tracer}(trs)
    diag_tr = merge(diag_tr, tr_default) # add PAR as default diagnostic

    ##### phytoplankton
    plank_name = keys(model.individuals.phytos)
    Nsp = length(plank_name)
    for i in 1:Nsp
        procs_sp = []
        for j in 1:nproc_phyto
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{phytoplankton}(procs_sp)
        procs_sp_d = []
        for j in 1:5
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp_d, proc)
        end
        diag_proc_default = NamedTuple{(:num, :graz, :mort, :dvid, :ptc)}(procs_sp_d)
        diag_proc = merge(diag_proc, diag_proc_default) # add num, graz, mort, and dvid as default diagnostics
        push!(phyto_procs, diag_proc)
    end
    diag_phyto = NamedTuple{plank_name}(phyto_procs)

    ##### abiotic particle
    abiotic_name = keys(model.individuals.abiotics)
    Nsa = length(abiotic_name)
    for i in 1:Nsa
        procs_sp = []
        for j in 1:nproc_abiotic
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{abiotic_particle}(procs_sp)
        procs_sp_d = []
        for j in 1:1
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp_d, proc)
        end
        diag_proc_default = NamedTuple{(:num,)}(procs_sp_d)
        diag_proc = merge(diag_proc, diag_proc_default) # add num as default diagnostics
        push!(abiotic_procs, diag_proc)
    end
    diag_abiotic = NamedTuple{abiotic_name}(abiotic_procs)

    ##### colony
    colony_name = keys(model.individuals.colonies)
    Ncl = length(colony_name)
    for i in 1:Ncl
        sp_name = keys(model.individuals.colonies[i].spcs)
        Nsp_cl = length(sp_name)
        sp_procs = []
        for j in 1:Nsp_cl
            procs_sp = []
            for k in 1:nproc_colony
                proc = zeros(FT, total_size) |> array_type(model.arch)
                push!(procs_sp, proc)
            end
            diag_proc = NamedTuple{colony}(procs_sp)
            procs_sp_d = []
            for k in 1:4
                proc = zeros(FT, total_size) |> array_type(model.arch)
                push!(procs_sp_d, proc)
            end
            diag_proc_default = NamedTuple{(:num, :graz, :mort, :dvid)}(procs_sp_d)
            diag_proc = merge(diag_proc, diag_proc_default) # add num, graz, mort, and dvid as default diagnostics
            push!(sp_procs, diag_proc)
        end
        diag_sp = NamedTuple{sp_name}(sp_procs)
        push!(colony_procs, diag_sp)
    end
    diag_colony = NamedTuple{colony_name}(colony_procs)
    
    diagnostics = PlanktonDiagnostics(diag_phyto, diag_abiotic, diag_colony, diag_tr, iteration_interval)

    return diagnostics
end

function show(io::IO, diags::PlanktonDiagnostics)
    if diags.abiotics == NamedTuple(;)
        s = "├── No abiotic particles available\n"
    else
        s = "├── diagnostics of abiotic particles: $(keys(diags.abiotics.sa1))\n"
    end
    if diags.colonies == NamedTuple(;)
        t = "├── No colony available\n"
    else
        t = "├── diagnostics of colony: $(keys(diags.colonies.cl1.sp1))\n"
    end
    print(io, "PlanktonDiagnostics:\n",
              "├── diagnostics of tracers: $(keys(diags.tracer))\n",
              "├── diagnostics of phytoplankton: $(keys(diags.phytos.sp1))\n",
              s,
              t,
              "└── save averaged diagnostics every $(diags.iteration_interval) timesteps")
end

function diag_avail(tracer, plank, abiotic, colony, model)
    tracer_avail = tracer_avail_diags()
    plank_avail  = (keys(components(model.individuals.phytos.sp1.data))..., :num)
    if keys(model.individuals.abiotics) == ()
        abiotic_avail = ()
    else
        abiotic_avail = (keys(components(model.individuals.abiotics.sa1.data))..., :num)
    end
    if keys(model.individuals.colonies) == ()
        colony_avail = ()
    else
        colony_avail = (keys(components(model.individuals.colonies.cl1.spcs.sp1.data))..., :num)
    end

    for i in eachindex(tracer)
        if tracer[i] ∉ tracer_avail
            throw(ArgumentError("$(tracer[i]) is not one of the diagnostics"))
        end
    end

    for i in eachindex(plank)
        if plank[i] ∉ plank_avail
            throw(ArgumentError("$(plank[i]) is not one of the diagnostics"))
        end
    end

    for i in eachindex(abiotic)
        if abiotic[i] ∉ abiotic_avail
            throw(ArgumentError("$(abiotic[i]) is not one of the diagnostics"))
        end
    end

    for i in eachindex(colony)
        if colony[i] ∉ colony_avail
            throw(ArgumentError("$(colony[i]) is not one of the diagnostics"))
        end
    end
end

function tracer_avail_diags()
    return (:PAR, :DIC, :DOC, :POC, :NH4, :NO3, :O2, :DON, :PON, :PO4, :DOP, :POP, :DFe, :PFe_inorg, :PFe_bio, :Dust)
end