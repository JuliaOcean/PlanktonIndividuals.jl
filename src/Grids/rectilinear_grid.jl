struct RectilinearGrid{FT, TX, TY, TZ, V, A3} <: AbstractGrid{FT, TX, TY, TZ}
    # corrdinates at cell centers, unit: meter
    xC::V
    yC::V
    zC::V
    # corrdinates at cell faces, unit: meter
    xF::V
    yF::V
    zF::V
    # grid spacing, unit: meter
    Δx::FT
    Δy::FT
    dzC::V
    dzF::V
    # number of grid points
    Nx::Int
    Ny::Int
    Nz::Int
    # number of halo points
    Hx::Int
    Hy::Int
    Hz::Int
    # landmask to indicate where is the land
    landmask::A3
end

"""
    RectilinearGrid(;size, x, y, z,
                     FT = Float32,
                     topology = (Periodic, Periodic, Bounded),
                     landmask = nothing,
                     halo = (2, 2, 2))
Creats a `RectilinearGrid` struct with `size = (Nx, Ny, Nz)` grid points.
    `x` and `y` directions must be regular spaced,
    `z` direction can be vertically stretched or regular spaced.

Keyword Arguments (Required)
============================
- `size` : A tuple prescribing the number of grid points. 
                `size` is a 3-tuple no matter for 3D, 2D, or 1D model.
- `x` and `y` : A 2-tuple that specify the start and end points of the domain.
- `z` : is either a (1) 1D array that specifies the locations of cell faces in z direction,
                or (2) 2-tuples that specify the start and end points of the domain.
                Vertical indexing starts from surface and use negative numbers for depth.

Keyword Arguments (Optional)
============================
- `FT`: Floating point data type. Default: `Float32`.
- `topology` : A 3-tuple specifying the topology of the domain.
                The topology can be either Periodic or Bounded in each direction.
- `landmask` : a 3-dimentional array to indicate where the land is.
- `halo` : A tuple of integers that specifies the size of the halo region of cells
                surrounding the physical interior for each direction.
                `halo` is a 3-tuple no matter for 3D, 2D, or 1D model.
                At least 2 halo points are needed for DST3FL advection scheme.
"""
function RectilinearGrid(;size, x, y, z,
                          FT = Float32,
                          topology = (Periodic, Periodic, Bounded),
                          landmask = nothing,
                          halo = (2, 2, 2))
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    TX, TY, TZ = validate_topology(topology)

    @assert isa(x, Tuple{<:Number, <:Number})
    @assert isa(y, Tuple{<:Number, <:Number})
    x₁, x₂ = x
    y₁, y₂ = y
    @assert x₁ < x₂
    @assert y₁ < y₂
    Δx = (x₂ - x₁) / Nx
    Δy = (y₂ - y₁) / Ny

    xF = collect(FT, range(-Hx * Δx, (Nx + Hx - 1) * Δx, length = Nx + 2 * Hx))
    yF = collect(FT, range(-Hy * Δy, (Ny + Hy - 1) * Δy, length = Ny + 2 * Hy))

    xC = collect(FT, range((0.5 - Hx) * Δx, (Nx + Hx - 0.5) * Δx, length = Nx + 2 * Hx))
    yC = collect(FT, range((0.5 - Hy) * Δy, (Ny + Hy - 0.5) * Δy, length = Ny + 2 * Hy))

    if isa(z, Tuple{<:Number, <:Number})
        z₁, z₂ = z
        z = collect(FT, range(z₁, z₂, length = Nz+1))
    elseif isa(z, AbstractVector)
        z₁, z₂ = z[1], z[end]
    else
        throw(ArgumentError("z must be a 2-tuple or 1D array"))
    end

    @assert z₁ > z₂
    @assert Base.length(z) == Nz + 1

    zF = zeros(FT, Nz+2Hz)
    zC = zeros(FT, Nz+2Hz)
    dzF = zeros(FT, Nz+2Hz)
    dzC = zeros(FT, Nz+2Hz)

    zF[1+Hz:Nz+Hz+1] .= FT.(z)
    zC[1+Hz:Nz+Hz] .= (zF[1+Hz:Nz+Hz] .+ zF[2+Hz:Nz+Hz+1]) ./ 2.0f0
    dzF[1+Hz:Nz+Hz] .= zF[1+Hz:Nz+Hz] .- zF[2+Hz:Nz+Hz+1]
    dzC[1+Hz:Nz+Hz-1] .= zC[1+Hz:Nz+Hz-1] .- zC[2+Hz:Nz+Hz]

    ##### fill halos
    @views @. dzF[1:Hz] = dzF[Hz+1]
    @views @. dzC[1:Hz] = dzF[Hz+1]
    @views @. dzF[Nz+Hz+1:Nz+2*Hz] = dzF[Nz+Hz]
    @views @. dzC[Nz+Hz:Nz+2*Hz] = dzF[Nz+Hz]

    @views @. zF[1:Hz] = zF[1+Hz:2*Hz] + dzF[1]*Hz
    @views @. zF[Nz+Hz+2:Nz+2*Hz] = zF[Nz+2:Nz+Hz] - dzF[end]*Hz
    for i in 1:Hz
        zC[Nz+Hz+i] = zC[Nz+Hz+i-1] - dzF[end]
    end
    for i in Hz:-1:1
        zC[i] = zC[i+1] + dzF[1]
    end

    landmask = landmask_validation(landmask, Nx, Ny, Nz, Hx, Hy, Hz, FT, TX, TY)

    return RectilinearGrid{FT, TX, TY, TZ, typeof(xF), typeof(landmask)}(
        xC, yC, zC, xF, yF, zF, Δx, Δy, dzC, dzF, Nx, Ny, Nz, Hx, Hy, Hz, landmask)
end

function show(io::IO, g::RectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    xL, xR = g.xF[g.Hx+1], g.xF[g.Hx+1+g.Nx]
    yL, yR = g.yF[g.Hy+1], g.yF[g.Hy+1+g.Ny]
    zL, zR = g.zF[g.Hz+1], g.zF[g.Hz+1+g.Nz]
    dzF_min = minimum(g.dzF)
    dzF_max = maximum(g.dzF)
    print(io, "RegularRectilinearGrid{$FT, $TX, $TY, $TZ}\n",
              "domain: x ∈ [$xL, $xR], y ∈ [$yL, $yR], z ∈ [$zL, $zR]\n",
              "topology (Tx, Ty, Tz):     ", (TX, TY, TZ), '\n',
              "resolution (Nx, Ny, Nz):   ", (g.Nx, g.Ny, g.Nz), '\n',
              "halo size (Hx, Hy, Hz):    ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): ", g.Δx, ", ", g.Δy, ", [min=", dzF_min, ", max=", dzF_max,"])")
end

function short_show(grid::RectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    return "RegularRectilinearGrid{$FT, $TX, $TY, $TZ}(Nx=$(grid.Nx), Ny=$(grid.Ny), Nz=$(grid.Nz))"
end


function validate_topology(topology)
    for T in topology
        if !isa(T(), AbstractTopology)
            e = "$T is not a valid topology!" * 
                "Valid topologies are: Periodic and Bounded."
            throw(ArgumentError(e))
        end
    end

        return topology
end