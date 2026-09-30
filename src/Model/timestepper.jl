function timestepper(arch::Architecture, FT::DataType, g::AbstractGrid, maxN, intac::Union{Nothing, AbstractArray}, palat::Palat)
    vel₀ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel½ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel₁ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))

    Gcs = tracers_init(arch, g, bgc_tracer_names, FT)
    tracer_temp = tracers_init(arch, g, bgc_tracer_names, FT)
    plk = tracers_init(arch, g, bgc_tracer_names, FT)

    par = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    par₀= zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    Chl = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    pop = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    flux_sink = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    temp = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    PARF = zeros(FT, g.Nx, g.Ny) |> array_type(arch)

    rnd = StructArray(x = zeros(FT, maxN), y = zeros(FT, maxN), z = zeros(FT, maxN))
    rnd_d = replace_storage(array_type(arch), rnd)

    rnd_3d = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)

    velos = StructArray(x  = zeros(FT, maxN), y  = zeros(FT, maxN), z  = zeros(FT, maxN),
                        u1 = zeros(FT, maxN), v1 = zeros(FT, maxN), w1 = zeros(FT, maxN),
                        u2 = zeros(FT, maxN), v2 = zeros(FT, maxN), w2 = zeros(FT, maxN),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    trs = StructArray(NH4 = zeros(FT, maxN), NO3 = zeros(FT, maxN), PO4 = zeros(FT, maxN), 
                      DOC = zeros(FT, maxN), DFe = zeros(FT, maxN), O2  = zeros(FT, maxN),
                      par = zeros(FT, maxN), 
                      T   = zeros(FT, maxN), pop = zeros(FT, maxN), dpar= zeros(FT, maxN), 
                      idc = zeros(FT, maxN), idc_int = zeros(Int, maxN))
    trs_d = replace_storage(array_type(arch), trs)

    ts = timestepper(Gcs, tracer_temp, vel₀, vel½, vel₁, PARF, temp, flux_sink, plk, par, par₀, Chl, pop, 
    rnd_d, rnd_3d, velos_d, trs_d, intac, palat)

    return ts
end
