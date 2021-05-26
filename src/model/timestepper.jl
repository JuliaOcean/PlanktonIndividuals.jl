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
    chl::AbstractArray  # a (Cu)Array to store Chl field of each timestep
    pop::AbstractArray  # a (Cu)Array to store population field of each timestep
    rnd::AbstractArray  # a StructArray of random numbers for plankton diffusion or grazing, mortality and division.
    velos::AbstractArray# a StructArray of intermediate values for RK4 particle advection
    nuts::AbstractArray # a StructArray of nutrients of each individual
end

function timestepper(arch::Architecture, g::AbstractGrid, N, maxN)
    vel₀ = (u = Field(arch, g), v = Field(arch, g), w = Field(arch, g))
    vel½ = (u = Field(arch, g), v = Field(arch, g), w = Field(arch, g))
    vel₁ = (u = Field(arch, g), v = Field(arch, g), w = Field(arch, g))

    Gcs = nutrients_init(arch, g)
    nut_temp = nutrients_init(arch, g)
    plk = nutrients_init(arch, g)

    par = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    chl = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    pop = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)

    temp = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    PARF = zeros(g.Nx, g.Ny) |> array_type(arch)

    rnd = StructArray(x = zeros(maxN), y = zeros(maxN), z = zeros(maxN))
    rnd_d = replace_storage(array_type(arch), rnd)

    velos = StructArray(x  = zeros(maxN), y  = zeros(maxN), z  = zeros(maxN),
                        u1 = zeros(maxN), v1 = zeros(maxN), w1 = zeros(maxN),
                        u2 = zeros(maxN), v2 = zeros(maxN), w2 = zeros(maxN),
                        u3 = zeros(maxN), v3 = zeros(maxN), w3 = zeros(maxN),
                        u4 = zeros(maxN), v4 = zeros(maxN), w4 = zeros(maxN),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    nuts = StructArray(NH4 = zeros(maxN), NO3 = zeros(maxN), PO4 = zeros(maxN), DOC = zeros(maxN),
                       par = zeros(maxN), T   = zeros(maxN), pop = zeros(maxN))
    nuts_d = replace_storage(array_type(arch), nuts)

    ts = timestepper(Gcs, nut_temp, vel₀, vel½, vel₁, PARF, temp, plk, par, chl, pop, rnd_d, velos_d, nuts_d)

    return ts
end