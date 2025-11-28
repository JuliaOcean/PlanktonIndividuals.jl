mutable struct timestepper
    Gcs::NamedTuple     # a NamedTuple same as tracers to store tendencies
    tracer_temp::NamedTuple# a NamedTuple same as tracers to store tracers fields in multi-dims advection scheme
    vel₀::NamedTuple    # a NamedTuple with u, v, w velocities
    vel½::NamedTuple    # a NamedTuple with u, v, w velocities
    vel₁::NamedTuple    # a NamedTuple with u, v, w velocities
    PARF::AbstractArray # a (Cu)Array to store surface PAR field of each timestep 
    temp::AbstractArray # a (Cu)Array to store temperature field of each timestep
    flux_sink::AbstractArray # a (Cu)Array to store sinking flux field of each timestep
    plk::NamedTuple     # a NamedTuple same as tracers to store interactions with individuals
    par::AbstractArray  # a (Cu)Array to store PAR field of each timestep
    par₀::AbstractArray # a (Cu)Array to store PAR field of the previous timestep
    Chl::AbstractArray  # a (Cu)Array to store Chl field of each timestep
    pop::AbstractArray  # a (Cu)Array to store population field of each timestep
    rnd::AbstractArray  # a StructArray of random numbers for plankton diffusion or grazing, mortality and division.
    rnd_3d::AbstractArray#a (Cu)Array of random numbers for tracer-particle interatcion
    velos::AbstractArray# a StructArray of intermediate values for RK4 particle advection
    trs::AbstractArray  # a StructArray of tracers of each individual
    intac::Union{Nothing, AbstractArray}   # Top-K candidate phyto IDs for abiotic particles
    palat::Palat        # a `Palat` to store the interaction between species
end

function timestepper(arch::Architecture, FT::DataType, g::AbstractGrid, maxN1, maxN2, intac::Union{Nothing,AbstractArray}, palat::Palat)
    vel₀ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel½ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel₁ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))

    Gcs = tracers_init(arch, g, FT)
    tracer_temp = tracers_init(arch, g, FT)
    plk = tracers_init(arch, g, FT)

    par = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    par₀= zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    Chl = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    pop = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    flux_sink = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    temp = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    PARF = zeros(FT, g.Nx, g.Ny) |> array_type(arch)

    maxN = max(maxN1, maxN2) 

    rnd = StructArray(x = zeros(FT, maxN), y = zeros(FT, maxN), z = zeros(FT, maxN))
    rnd_d = replace_storage(array_type(arch), rnd)

    rnd_3d = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)

    velos = StructArray(x  = zeros(FT, maxN), y  = zeros(FT, maxN), z  = zeros(FT, maxN),
                        u1 = zeros(FT, maxN), v1 = zeros(FT, maxN), w1 = zeros(FT, maxN),
                        u2 = zeros(FT, maxN), v2 = zeros(FT, maxN), w2 = zeros(FT, maxN),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    trs = StructArray(NH4 = zeros(FT, maxN), NO3 = zeros(FT, maxN), PO4 = zeros(FT, maxN), 
                      DOC = zeros(FT, maxN), DFe = zeros(FT, maxN), par = zeros(FT, maxN), 
                      T   = zeros(FT, maxN), pop = zeros(FT, maxN), dpar= zeros(FT, maxN), 
                      idc = zeros(FT, maxN), idc_int = zeros(Int, maxN))
    trs_d = replace_storage(array_type(arch), trs)

    
    ts = timestepper(Gcs, tracer_temp, vel₀, vel½, vel₁, PARF, temp, flux_sink, plk, par, par₀, Chl, pop, rnd_d, rnd_3d, velos_d, trs_d, intac, palat)

    return ts
end