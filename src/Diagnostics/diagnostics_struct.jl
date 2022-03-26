mutable struct PlanktonDiagnostics
    plankton::NamedTuple       # for each species
    tracer::NamedTuple         # for tracers
    time_interval::Int64       # time interval that the diagnostics is time averaged
end

"""
    PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                        plankton=(:num, :graz, :mort, :dvid),
                        time_interval = 3600)

Generate a `PlanktonDiagnostics` structure.

Keyword Arguments (Optional)
============================
- `tracer` : a `Tuple` containing the names of nutrient fields to be diagnosed.
- `plankton` : a `Tuple` containing the names of physiological processes of plankton individuals to be diagnosed.
- `time_interval` : The time interval that diagnostics is averaged, an hour (3600 seconds) by default.
"""
function PlanktonDiagnostics(model; tracer=(:PAR,),
                            plankton=(:num, :graz, :mort, :dvid),
                            time_interval::Int64 = 3600)
    
    @assert isa(tracer, Tuple)
    @assert isa(plankton, Tuple)

    diag_avail(tracer, plankton, model.mode)

    ntr   = length(tracer)
    nproc = length(plankton)
    trs   = []
    procs = []

    total_size = (model.grid.Nx+model.grid.Hx*2, model.grid.Ny+model.grid.Hy*2, model.grid.Nz+model.grid.Hz*2)

    for i in 1:ntr
        tr = zeros(total_size) |> array_type(model.arch)
        push!(trs, tr)
    end
    tr_d = zeros(total_size) |> array_type(model.arch)
    tr_default = (PAR = tr_d,)

    diag_tr = NamedTuple{tracer}(trs)
    diag_tr = merge(diag_tr, tr_default) # add PAR as default diagnostic

    plank_name = keys(model.individuals.phytos)
    Nsp = length(plank_name)

    for j in 1:Nsp
        procs_sp = []
        for k in 1:nproc
            proc = zeros(total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{plankton}(procs_sp)

        procs_sp_d = []
        for l in 1:4
            proc = zeros(total_size) |> array_type(model.arch)
            push!(procs_sp_d, proc)
        end
        diag_proc_default = NamedTuple{(:num, :graz, :mort, :dvid)}(procs_sp_d)

        diag_proc = merge(diag_proc, diag_proc_default) # add num, graz, mort, and dvid as default diagnostics

        push!(procs, diag_proc)
    end
    diag_sp = NamedTuple{plank_name}(procs)

    diagnostics = PlanktonDiagnostics(diag_sp, diag_tr, time_interval)

    return diagnostics
end

function show(io::IO, diags::PlanktonDiagnostics)
    print(io, "PlanktonDiagnostics:\n",
              "├── diagnostics of tracers: $(keys(diags.tracer))\n",
              "├── diagnostics of individuals: $(keys(diags.plankton.sp1))\n",
              "└── save averaged diagnostics every $(diags.time_interval) seconds")
end

function diag_avail(tracer, plank, mode::AbstractMode)
    tracer_avail = tracer_avail_diags()
    plank_avail  = plank_avail_diags(mode)
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
end

function tracer_avail_diags()
    return (:PAR, :DIC, :DOC, :POC, :NH4, :NO3, :DON, :PON, :PO4, :DOP, :POP)
end

function plank_avail_diags(mode::AbstractMode)
    plank_avail = (:num, :graz, :mort, :dvid, :PS, :resp, :Bm, :Chl, :Th)
    if isa(mode, CarbonMode)
        nothing
    elseif isa(mode, QuotaMode)
        plank_avail = (:num, :graz, :mort, :dvid, :PS, :BS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :exu, :Bm, :Cq, :Nq, :Pq, :Chl)
    elseif isa(mode, MacroMolecularMode)
        plank_avail = (:num, :graz, :mort, :dvid, :PS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :ρChl, :S_PRO, :S_DNA, :S_RNA, :exu, :CH, :NST, :PST, :PRO, :DNA, :RNA, :Chl)
    end
    return plank_avail
end
