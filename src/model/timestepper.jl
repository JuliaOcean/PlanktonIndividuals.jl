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

    rnd = StructArray(x = zeros(λ*N), y = zeros(λ*N), z = zeros(λ*N))
    rnd_d = replace_storage(array_type(arch), rnd)

    velos = StructArray(x  = zeros(λ*N), y  = zeros(λ*N), z  = zeros(λ*N),
                        u1 = zeros(λ*N), v1 = zeros(λ*N), w1 = zeros(λ*N),
                        u2 = zeros(λ*N), v2 = zeros(λ*N), w2 = zeros(λ*N),
                        u3 = zeros(λ*N), v3 = zeros(λ*N), w3 = zeros(λ*N),
                        u4 = zeros(λ*N), v4 = zeros(λ*N), w4 = zeros(λ*N),
                        )
    velos_d = replace_storage(array_type(arch), velos)

    nuts = StructArray(NH4 = zeros(λ*N), NO3 = zeros(λ*N), PO4 = zeros(λ*N), DOC = zeros(λ*N),
                       αI  = zeros(λ*N), Tem = zeros(λ*N), pop = zeros(λ*N))
    nuts_d = replace_storage(array_type(arch), nuts)

    ts = timestepper(Gcs, MD1, MD2, MD3, vel₀, vel½, vel₁, plk, par, chl, pop, rnd_d, velos_d, nuts_d)

    return ts
end