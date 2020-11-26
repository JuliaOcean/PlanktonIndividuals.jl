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
    tmp::AbstractArray  # a (Cu)Array to store temporary individuals of each timestep
    rnd::AbstractArray  # a StructArray of random numbers for plankton diffusion or grazing, mortality and division.
    velos::AbstractArray# a StructArray of intermediate values for RK4 particle advection
    nuts::AbstractArray # a StructArray of nutrients of each individual
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

    par = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    chl = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    pop = zeros(g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)

    tmp = StructArray(x   = zeros(4N), y   = zeros(4N), z   = zeros(4N),
                      xi  = zeros(4N), yi  = zeros(4N), zi  = zeros(4N), 
                      iS  = zeros(4N), Sz  = zeros(4N), Bm  = zeros(4N), 
                      Cq  = zeros(4N), Nq  = zeros(4N), Pq  = zeros(4N), 
                      chl = zeros(4N), gen = zeros(4N), age = zeros(4N), 
                      ac  = zeros(4N), idx = zeros(4N),
                      graz= zeros(4N), mort= zeros(4N), dvid= zeros(4N))
    tmp_d = replace_storage(array_type(arch), tmp)

    rnd = StructArray(x = zeros(4N), y = zeros(4N), z = zeros(4N))
    rnd_d = replace_storage(array_type(arch), rnd)

    velos = StructArray(x  = zeros(4N), y  = zeros(4N), z  = zeros(4N),
                        u1 = zeros(4N), v1 = zeros(4N), w1 = zeros(4N),
                        u2 = zeros(4N), v2 = zeros(4N), w2 = zeros(4N),
                        u3 = zeros(4N), v3 = zeros(4N), w3 = zeros(4N),
                        u4 = zeros(4N), v4 = zeros(4N), w4 = zeros(4N),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    nuts = StructArray(NH4 = zeros(4N), NO3 = zeros(4N), PO4 = zeros(4N), DOC = zeros(4N),
                       αI  = zeros(4N), Tem = zeros(4N), pop = zeros(4N))
    nuts_d = replace_storage(array_type(arch), nuts)

    ts = timestepper(Gcs, MD1, MD2, MD3, vel₀, vel½, vel₁, plk,
                     par, chl, pop, tmp_d, rnd_d, velos_d, nuts_d)

    return ts
end