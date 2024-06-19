mutable struct timestepper
    Gcs::NamedTuple     # a NamedTuple same as nutrients to store tendencies
    nut_temp::NamedTuple# a NamedTuple same as nutrients to store nutrients fields in multi-dims advection scheme
    vel₀::NamedTuple    # a NamedTuple with u, v, w velocities
    vel½::NamedTuple    # a NamedTuple with u, v, w velocities
    vel₁::NamedTuple    # a NamedTuple with u, v, w velocities
    PARF::AbstractArray # a (Cu)Array to store surface PAR field of each timestep 
    temp::AbstractArray # a (Cu)Array to store temperature field of each timestep
    plk::NamedTuple     # a NamedTuple same as nutrients to store interactions with individuals
    par::AbstractArray  # a (Cu)Array to store PAR field of each timestep
    Chl::AbstractArray  # a (Cu)Array to store Chl field of each timestep
    pop::AbstractArray  # a (Cu)Array to store population field of each timestep
    rnd::AbstractArray  # a StructArray of random numbers for plankton diffusion or grazing, mortality and division.
    velos::AbstractArray# a StructArray of intermediate values for RK4 particle advection
    nuts::AbstractArray # a StructArray of nutrients of each individual
end

function timestepper(arch::Architecture, FT::DataType, g::AbstractGrid, maxN)
    vel₀ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel½ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel₁ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))

    Gcs = nutrients_init(arch, g, FT)
    nut_temp = nutrients_init(arch, g, FT)
    plk = nutrients_init(arch, g, FT)

    par = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    Chl = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    pop = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)

    temp = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    PARF = zeros(FT, g.Nx, g.Ny) |> array_type(arch)

    rnd = StructArray(x = zeros(FT, maxN), y = zeros(FT, maxN), z = zeros(FT, maxN))
    rnd_d = replace_storage(array_type(arch), rnd)

    velos = StructArray(x  = zeros(FT, maxN), y  = zeros(FT, maxN), z  = zeros(FT, maxN),
                        u1 = zeros(FT, maxN), v1 = zeros(FT, maxN), w1 = zeros(FT, maxN),
                        u2 = zeros(FT, maxN), v2 = zeros(FT, maxN), w2 = zeros(FT, maxN),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    nuts = StructArray(NH4 = zeros(FT, maxN), NO3 = zeros(FT, maxN), PO4 = zeros(FT, maxN), 
                       DOC = zeros(FT, maxN), FeT = zeros(FT, maxN), idc = zeros(Int,maxN),
                       par = zeros(FT, maxN), T   = zeros(FT, maxN), pop = zeros(FT, maxN))
    nuts_d = replace_storage(array_type(arch), nuts)

    ts = timestepper(Gcs, nut_temp, vel₀, vel½, vel₁, PARF, temp, plk, par, Chl, pop, rnd_d, velos_d, nuts_d)

    return ts
end