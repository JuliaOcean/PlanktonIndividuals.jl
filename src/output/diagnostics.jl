mutable struct PlanktonDiagnostics
    plankton::NamedTuple       # for each species
    tracer::NamedTuple         # for tracers
    frequency::Int64           # frequency of diagnostics (in number of time steps)
end

"""
    PlanktonDiagnostics(model;
                        tracer=(:PAR, :NH4, :NO3, :DOC),
                        plankton=(:num, :graz, :mort, :dvid),
                        frequency = 1)

Generate a `PlanktonDiagnostics` structure.

Keyword Arguments (Optional)
============================
- `tracer` : a `Tuple` containing the names of nutrient fields to be diagnosed.
- `plankton` : a `Tuple` containing the names of physiological processes of plankton individuals to be diagnosed.
- `frequency` : frequency of diagnostics (in numbers of time steps), diagnose every time step by default.
"""

function PlanktonDiagnostics(model;
                            tracer=(:PAR, :NH4, :NO3, :DOC),
                            plankton=(:num, :graz, :mort, :dvid),
                            frequency::Int64 = 1)
    ntr   = length(tracer)
    nproc = length(plankton)
    trs   = []
    procs = []

    total_size = (model.grid.Nx+model.grid.Hx*2, model.grid.Ny+model.grid.Hy*2, model.grid.Nz+model.grid.Hz*2)

    for i in 1:ntr
        tr = zeros(total_size) |> array_type(model.arch)
        push!(trs, tr)
    end

    diag_tr = NamedTuple{tracer}(trs)

    Nsp = length(model.individuals.phytos)

    for j in 1:Nsp
        procs_sp = []
        for k in 1:nproc
            proc = zeros(total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{plankton}(procs_sp)
        push!(procs, diag_proc)
    end
    plank_name = plank_names[1:Nsp]
    diag_sp = NamedTuple{plank_name}(procs)

    diagnostics = PlanktonDiagnostics(diag_sp, diag_tr, frequency)

    return diagnostics
end

function show(io::IO, diags::PlanktonDiagnostics)
    Nsp = length(model.individuals.phytos)
    N = Int(dot(model.individuals.phytos.sp1.data.ac,model.individuals.phytos.sp1.data.ac))
    cap = length(model.individuals.phytos.sp1.data.ac)
    print(io, "diagnostics of tracers: $(keys(diags.tr))\n",
              "diagnostics of individuals: $(keys(diags.spcs.sp1))",
              "save averaged diagnostics every $(diags.frequency) time steps\n")
end