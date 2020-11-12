mutable struct timestepper
    Gcs::NamedTuple     # a NamedTuple same as nutrients to store tendencies
    MD1::NamedTuple     # a NamedTuple same as nutrients to store nutrients fields in multi-dims advection scheme
    MD2::NamedTuple     # a NamedTuple same as nutrients to store nutrients fields in multi-dims advection scheme
    MD3::NamedTuple     # a NamedTuple same as nutrients to store nutrients fields in multi-dims advection scheme
    vel₀::NamedTuple    # a NamedTuple with u, v, w velocities
    vel½::NamedTuple    # a NamedTuple with u, v, w velocities
    vel₁::NamedTuple    # a NamedTuple with u, v, w velocities
    plk::NamedTuple     # a NamedTuple same as nutrients to store interactions with individuals
    par::AbstractArray  # a (Cu)Array to store PAR field of each timestep
    chl::AbstractArray  # a (Cu)Array to store Chl field of each timestep
    pop::AbstractArray  # a (Cu)Array to store population field of each timestep
    cts::NamedTuple     # a NamedTuple to store counts of each timestep
    tmp::AbstractArray  # a (Cu)Array to store temporary individuals of each timestep
    rnd::AbstractArray  # a StructArray of random numbers for plankton diffusion or grazing, mortality and division.
    velos::AbstractArray# a StructArray of intermediate values for RK4 particle advection
    coord::AbstractArray# a StructArray of coordinates of each individual
    nuts::AbstractArray # a StructArray of nutrients of each individual
    proc::AbstractArray # a StructArray of process fluxes of each individual
end

function timestepper(arch::Architecture, g::Grids, N)
    vel₀ = (u = Field(arch, g), v = Field(arch, g), w = Field(arch, g))
    vel½ = (u = Field(arch, g), v = Field(arch, g), w = Field(arch, g))
    vel₁ = (u = Field(arch, g), v = Field(arch, g), w = Field(arch, g))

    Gcs = nutrients_init(arch, g)
    MD1 = nutrients_init(arch, g)
    MD2 = nutrients_init(arch, g)
    MD3 = nutrients_init(arch, g)
    plk = nutrients_init(arch, g)

    par = zeros(g.Nx, g.Ny, g.Nz) |> array_type(arch)
    chl = zeros(g.Nx, g.Ny, g.Nz) |> array_type(arch)
    pop = zeros(g.Nx, g.Ny, g.Nz) |> array_type(arch)

    cts = (chl = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           pop = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           DIC = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           NH4 = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           NO3 = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           PO4 = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           DOC = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           POC = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           DON = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           PON = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           DOP = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch),
           POP = zeros(g.Nx, g.Ny, g.Nz, 4N) |> array_type(arch))

    tmp = zeros(4N,60) |> array_type(arch)

    rnd = StructArray(x = zeros(4N), y = zeros(4N), z = zeros(4N))
    rnd_d = replace_storage(array_type(arch), rnd)

    velos = StructArray(x = zeros(4N), y = zeros(4N), z = zeros(4N),
                        u⁺= zeros(4N), v⁺= zeros(4N), w⁺= zeros(4N),
                        u⁻= zeros(4N), v⁻= zeros(4N), w⁻= zeros(4N),
                        u1= zeros(4N), v1= zeros(4N), w1= zeros(4N),
                        u2= zeros(4N), v2= zeros(4N), w2= zeros(4N),
                        u3= zeros(4N), v3= zeros(4N), w3= zeros(4N),
                        u4= zeros(4N), v4= zeros(4N), w4= zeros(4N),
                        xd= zeros(4N), yd= zeros(4N), zd= zeros(4N),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    coord = StructArray(x = zeros(4N), y = zeros(4N), z = zeros(4N))
    coord_d = replace_storage(array_type(arch), coord)

    nuts = StructArray(NH4 = zeros(4N), NO3 = zeros(4N), PO4 = zeros(4N), DOC = zeros(4N),
                       αI  = zeros(4N), Tem = zeros(4N), pop = zeros(4N))
    nuts_d = replace_storage(array_type(arch), nuts)

    proc = StructArray(PS   = zeros(4N), VDOC = zeros(4N), VNH4 = zeros(4N), VNO3 = zeros(4N),
                       VPO4 = zeros(4N), ρchl = zeros(4N), resp = zeros(4N), BS   = zeros(4N),
                       exu  = zeros(4N), grz  = zeros(4N), mort = zeros(4N), dvid = zeros(4N))
    proc_d = replace_storage(array_type(arch), proc)

    ts = timestepper(Gcs, MD1, MD2, MD3, vel₀, vel½, vel₁, plk,
                     par, chl, pop, cts, tmp, rnd_d, velos_d, coord_d, nuts_d, proc_d)

    return ts
end