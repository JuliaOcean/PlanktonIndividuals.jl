mutable struct PlanktonDiagnostics
    phytoplankton::NamedTuple       # for each species of phytoplankton
    abiotic_particle::NamedTuple    # for each species of abiotic particle
    tracer::NamedTuple              # for tracers
    iteration_interval::Int         # time interval that the diagnostics is time averaged
end

"""
    PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                               phytoplankton = (:num, :graz, :mort, :dvid),
                               abiotic_particle = (:num),
                               time_interval = 1)

Generate a `PlanktonDiagnostics` structure.

Keyword Arguments (Optional)
============================
- `tracer` : a `Tuple` containing the names of nutrient fields to be diagnosed.
- `phytoplankton` : a `Tuple` containing the names of physiological processes of phytoplankton individuals to be diagnosed.
- `abiotic_particle` : a `Tuple` containing the names of state variables of abiotic particles to be diagnosed.
- `iteration_interval` : The number of timesteps that diagnostics is averaged, 1 iteration by default.
"""
function PlanktonDiagnostics(model; tracer=(),
                                    phytoplankton = (:num, :graz, :mort, :dvid),
                                    abiotic_particle = (:num,),
                                    iteration_interval::Int = 1)
    
    @assert isa(tracer, Tuple)
    @assert isa(phytoplankton, Tuple)
    @assert isa(abiotic_particle, Tuple)

    diag_avail(tracer, phytoplankton, abiotic_particle, model.mode)

    ntr   = length(tracer)
    nproc_phyto = length(phytoplankton)
    nproc_abiotic = length(abiotic_particle)
    trs   = []
    phyto_procs = []
    abiotic_procs = []
    FT = model.FT

    total_size = (model.grid.Nx+model.grid.Hx*2, model.grid.Ny+model.grid.Hy*2, model.grid.Nz+model.grid.Hz*2)

    for i in 1:ntr
        tr = zeros(FT, total_size) |> array_type(model.arch)
        push!(trs, tr)
    end
    tr_d1 = zeros(FT, total_size) |> array_type(model.arch)
    tr_d2 = zeros(FT, total_size) |> array_type(model.arch)
    tr_default = (PAR = tr_d1, T = tr_d2)

    diag_tr = NamedTuple{tracer}(trs)
    diag_tr = merge(diag_tr, tr_default) # add PAR as default diagnostic

    plank_name = keys(model.individuals.phytos)
    Nsp_phyto = length(plank_name)

    for j in 1:Nsp_phyto
        procs_sp = []
        for k in 1:nproc_phyto
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{phytoplankton}(procs_sp)

        procs_sp_d = []
        for l in 1:4
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp_d, proc)
        end
        diag_proc_default = NamedTuple{(:num, :graz, :mort, :dvid)}(procs_sp_d)

        diag_proc = merge(diag_proc, diag_proc_default) # add num, graz, mort, and dvid as default diagnostics

        push!(phyto_procs, diag_proc)
    end
    diag_phyto = NamedTuple{plank_name}(phyto_procs)

    abiotic_name = keys(model.individuals.abiotics)
    Nsp_abiotic = length(abiotic_name)

    for k in 1:Nsp_abiotic
        procs_sp = []
        for l in 1:nproc_abiotic
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{abiotic_particle}(procs_sp)

        procs_sp_d = []
        for m in 1:1
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp_d, proc)
        end
        diag_proc_default = NamedTuple{(:num,)}(procs_sp_d)

        diag_proc = merge(diag_proc, diag_proc_default) # add num as default diagnostics

        push!(abiotic_procs, diag_proc)
    end
    diag_abiotic = NamedTuple{abiotic_name}(abiotic_procs)

    diagnostics = PlanktonDiagnostics(diag_phyto, diag_abiotic, diag_tr, iteration_interval)

    return diagnostics
end

function show(io::IO, diags::PlanktonDiagnostics)
    print(io, "PlanktonDiagnostics:\n",
              "├── diagnostics of tracers: $(keys(diags.tracer))\n",
              "├── diagnostics of phytoplankton: $(keys(diags.phytoplankton.sp1))\n",
              "├── diagnostics of abiotic particles: $(keys(diags.abiotic_particle.sp1))\n",
              "└── save averaged diagnostics every $(diags.iteration_interval) timesteps")
end

function diag_avail(tracer, plank, abiotic, mode::AbstractMode)
    tracer_avail = tracer_avail_diags()
    plank_avail  = plank_avail_diags(mode)
    abiotic_avail = abiotic_avail_diags()
    for i in 1:length(tracer)
        if length(findall(x->x==tracer[i], tracer_avail)) == 0
            throw(ArgumentError("$(tracer[i]) is not one of the diagnostics"))
        end
    end

    for i in 1:length(plank)
        if length(findall(x->x==plank[i], plank_avail)) == 0
            throw(ArgumentError("$(plank[i]) is not one of the diagnostics"))
        end
    end

    for i in 1:length(abiotic)
        if length(findall(x->x==abiotic[i], abiotic_avail)) == 0
            throw(ArgumentError("$(abiotic[i]) is not one of the diagnostics"))
        end
    end
end

function tracer_avail_diags()
    return (:PAR, :DIC, :DOC, :POC, :NH4, :NO3, :DON, :PON, :PO4, :DOP, :POP, :FeT, :DOFe, :POFe, :CHO)
end

function plank_avail_diags(mode::AbstractMode)
    plank_avail = (:num, :graz, :mort, :dvid, :PS, :Bm, :Chl)
    if isa(mode, CarbonMode)
        plank_avail = (:num, :graz, :mort, :dvid, :PS, :BS, :RP, :TD, :RS, :Bm, :Bd, :Chl)
    elseif isa(mode, QuotaMode)
        plank_avail = (:num, :graz, :mort, :dvid, :PS, :BS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :exu, :Bm, :Cq, :Nq, :Pq, :Chl)
    elseif isa(mode, MacroMolecularMode)
        plank_avail = (:num, :graz, :mort, :dvid, :PS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :ρChl, :S_PRO, :S_DNA, :S_RNA, :exu, :CH, :NST, :PST, :PRO, :DNA, :RNA, :Chl)
    elseif isa(mode, IronEnergyMode)
        plank_avail = (:num, :graz, :mort, :dvid, :PS, :CF, :ECF, :RS, :ERS, :NR, :ENR, :NF, :ENF, :BS, :VNH4, :VNO3, :VPO4, :VFe, :PS2ST, :ST2PS, :NR2ST, :ST2NR, :NF2ST, :ST2NF, :Bm, :exEn, :CH, :qNO3, :qNH4, :qP, :qFe, :qFePS, :qFeNR, :qFeNF, :Chl, :tdark)
    end
    return plank_avail
end

function abiotic_avail_diags()
    return (:num, :CHO, :CHOe, :ADS, :DEC)
end
