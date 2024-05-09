struct LatLonGrid{FT, TX, TY, TZ, V, A2, A3} <: AbstractGrid{FT, TX, TY, TZ}
    # corrdinates at cell centers, unit: degree
    xC::V
    yC::V
    # corrdinates at cell centers, unit: meter
    zC::V
    # corrdinates at cell faces, unit: degree
    xF::V
    yF::V
    # corrdinates at cell faces, unit: meter
    zF::V
    # grid spacing, unit: degree
    Δx::FT
    Δy::FT
    # grid spacing from center to center, unit: meter
    dxC::A2
    dyC::A2
    dzC::V
    # grid spacing from face to face, unit: meter
    dxF::A2
    dyF::A2
    dzF::V
    # areas and volume, unit: m² or m³
    Ax::A3
    Ay::A3
    Az::A2
    Vol::A3
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
    LatLonGrid(;size, lat, lon, z,
                FT = Float32,
                radius = 6370.0e3,
                landmask = nothing,
                halo = (2, 2, 2))
Creats a `LatLonGrid` struct with `size = (Nx, Ny, Nz)` grid points.

Keyword Arguments (Required)
============================
- `size` : A tuple prescribing the number of grid points. 
                `size` is a 3-tuple no matter for 3D, 2D, or 1D model.
- `lat` : A 2-tuple specifying the startind and ending points in latitudinal direction.
                Possible values are from -80 (80S) to 80 (80N).
- `lon` : A 2-tuple specifying the startind and ending points in longitudinal direction.
                Possible values are from -180 (180W) to 180 (180E).
- `z` : is either a (1) 1D array that specifies the locations of cell faces in z direction,
                or (2) 2-tuples that specify the start and end points of the domain.
                Vertical indexing starts from surface and use negative numbers for depth.

Keyword Arguments (Optional)
============================
- `FT`: Floating point data type. Default: `Float32`.
- `radius` : Specify the radius of the Earth used in the model, 6370.0e3 meters by default.
- `landmask` : a 3-dimentional array to indicate where the land is.
- `halo` : A tuple of integers that specifies the size of the halo region of cells
                surrounding the physical interior for each direction.
                `halo` is a 3-tuple no matter for 3D, 2D, or 1D model.
                At least 2 halo points are needed for DST3FL advection scheme.
"""
function LatLonGrid(;size, lat, lon, z,
                     FT = Float32,
                     radius = 6370.0e3,
                     landmask = nothing,
                     halo = (2, 2, 2))
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    lat₁, lat₂ = lat
    lon₁, lon₂ = lon

    if isa(z, Tuple{<:Number, <:Number})
        z₁, z₂ = z
        z = collect(FT, range(z₁, z₂, length = Nz+1))
    elseif isa(z, AbstractVector)
        z₁, z₂ = z[1], z[end]
    else
        throw(ArgumentError("z must be a 2-tuple or 1D array"))
    end

    @assert -180 <= lon₁ < lon₂ <= 180
    @assert -80 <= lat₁ < lat₂ <= 80
    @assert z₁ > z₂
    @assert Base.length(z) == Nz + 1

    TX = lon₁ == -180 && lon₂ == 180 ? Periodic : Bounded
    TY = Bounded
    TZ = Bounded

    Δx = (lon₂ - lon₁) / Nx
    Δy = (lat₂ - lat₁) / Ny

    xF = collect(FT, range(lon₁ - Hx * Δx, lon₁ + (Nx + Hx - 1) * Δx, length = Nx + 2 * Hx))
    yF = collect(FT, range(lat₁ - Hy * Δy, lat₁ + (Ny + Hy - 1) * Δy, length = Ny + 2 * Hy))

    xC = collect(FT, range(lon₁ + (0.5 - Hx) * Δx, lon₁ + (Nx + Hx - 0.5) * Δx, length = Nx + 2 * Hx))
    yC = collect(FT, range(lat₁ + (0.5 - Hy) * Δy, lat₁ + (Ny + Hy - 0.5) * Δy, length = Ny + 2 * Hy))

    # inclue halo points
    zF = zeros(FT, Nz+2Hz)
    zC = zeros(FT, Nz+2Hz)
    dzF = zeros(FT, Nz+2Hz)
    dzC = zeros(FT, Nz+2Hz)
    dxC = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    dyC = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    dxF = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    dyF = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    Ax = zeros(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    Ay = zeros(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    Az = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    Vol = zeros(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)

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

    for i in 1:Nx+2*Hx
        for j in 1:Ny+2*Hy
            dxC[i, j] = radius * cos(yC[j]*π/180) * deg2rad(Δx)
            dyC[i, j] = radius * deg2rad(Δy)
            dxF[i, j] = radius * cos(yC[j]*π/180) * deg2rad(Δx)
            dyF[i, j] = radius * deg2rad(Δy)
        end
    end
    for i in 1:Nx+2*Hx-1
        for j in 1:Ny+2*Hy-1
            Az[i,j]   = radius^2 * deg2rad(Δx) * (sin(yF[j+1]*π/180) - sin(yF[j]*π/180))
        end
    end
    Az[Nx+2*Hx, :] .= Az[Nx+2*Hx-1, :] 
    Az[1:Nx+2*Hx-1, Ny+2*Hy] .= Az[1:Nx+2*Hx-1, Ny+2*Hy-1] 

    for i in 1:Nx+2*Hx
        for j in 1:Ny+2*Hy
            for k in 1:Nz+2*Hz
                Ax[i,j,k]  = dzF[k] * dyF[i,j]
                Ay[i,j,k]  = dzF[k] * dxF[i,j]
                Vol[i,j,k] = Az[i,j] * dzF[k]
            end
        end
    end

    landmask = landmask_validation(landmask, Nx, Ny, Nz, Hx, Hy, Hz, FT, TX, TY)

    return LatLonGrid{FT, TX, TY, TZ, typeof(xF), typeof(dxC), typeof(Vol)}(
        xC, yC, zC, xF, yF, zF, Δx, Δy, dxC, dyC, dzC, dxF, dyF, dzF, Ax, Ay, Az, Vol, Nx, Ny, Nz, Hx, Hy, Hz, landmask)
end

"""
    LoadLatLonGrid(;grid_info, size, lat, lon,
                    FT = Float32,
                    landmask = nothing,
                    halo=(2,2,2))
Creats a `LatLonGrid` struct with `size = (Nx, Ny, Nz)` grid points.

Keyword Arguments (Required)
============================
- `grid_info` : A NamedTuple contains external grid information (e.g. from MITgcm), please
                    refer to documentation for the required format.
- `size` : A tuple prescribing the number of grid points. 
                `size` is a 3-tuple no matter for 3D, 2D, or 1D model.
- `lat` : A 2-tuple specifying the startind and ending points in latitudinal direction.
                Possible values are from -80 (80S) to 80 (80N).
- `lon` : A 2-tuple specifying the startind and ending points in longitudinal direction.
                Possible values are from -180 (180W) to 180 (180E).

Keyword Arguments (Optional)
============================
- `FT`: Floating point data type. Default: `Float32`.
- `landmask` : a 3-dimentional array to indicate where the land is.
- `halo` : A tuple of integers that specifies the size of the halo region of cells
                surrounding the physical interior for each direction.
                `halo` is a 3-tuple no matter for 3D, 2D, or 1D model.
                At least 2 halo points are needed for DST3FL advection scheme.
"""
function LoadLatLonGrid(;grid_info, size, lat, lon, 
                         FT = Float32, landmask = nothing, halo=(2,2,2))
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    lat₁, lat₂ = lat
    lon₁, lon₂ = lon

    @assert -180 <= lon₁ < lon₂ <= 180
    @assert -80 <= lat₁ < lat₂ <= 80

    TX = lon₁ == -180 && lon₂ == 180 ? Periodic : Bounded
    TY = Bounded
    TZ = Bounded

    Δx = (lon₂ - lon₁) / Nx
    Δy = (lat₂ - lat₁) / Ny

    xF = collect(FT, range(lon₁ - Hx * Δx, lon₁ + (Nx + Hx - 1) * Δx, length = Nx + 2 * Hx))
    yF = collect(FT, range(lat₁ - Hy * Δy, lat₁ + (Ny + Hy - 1) * Δy, length = Ny + 2 * Hy))

    xC = collect(FT, range(lon₁ + (0.5f0 - Hx) * Δx, lon₁ + (Nx + Hx - 0.5f0) * Δx, length = Nx + 2 * Hx))
    yC = collect(FT, range(lat₁ + (0.5f0 - Hy) * Δy, lat₁ + (Ny + Hy - 0.5f0) * Δy, length = Ny + 2 * Hy))

    zF = zeros(FT, Nz+2Hz)
    zC = zeros(FT, Nz+2Hz)
    dzF = zeros(FT, Nz+2Hz)
    dzC = zeros(FT, Nz+2Hz)
    dxC = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    dyC = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    dxF = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    dyF = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    Ax = zeros(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    Ay = zeros(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    Az = zeros(FT, Nx+2*Hx, Ny+2*Hy)
    Vol = zeros(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    hFW = ones(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    hFS = ones(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)
    hFC = ones(FT, Nx+2*Hx, Ny+2*Hy, Nz+2*Hz)

    zF[1+Hz:Nz+Hz+1] = FT.(grid_info.RF)
    zC[1+Hz:Nz+Hz] = FT.(grid_info.RC)
    dzF[1+Hz:Nz+Hz] = FT.(grid_info.DRF)
    dzC[1+Hz:Nz+Hz] = FT.(grid_info.DRC);  dzC[1+Hz] *= 2.0f0 # First laryer only has half of the grid
    dxC[1+Hx:Nx+Hx, 1+Hy:Ny+Hy] = FT.(grid_info.DXC)
    dxF[1+Hx:Nx+Hx, 1+Hy:Ny+Hy] = FT.(grid_info.DXG)
    dyC[1+Hx:Nx+Hx, 1+Hy:Ny+Hy] = FT.(grid_info.DYC)
    dyF[1+Hx:Nx+Hx, 1+Hy:Ny+Hy] = FT.(grid_info.DYG)
    Az[1+Hx:Nx+Hx, 1+Hy:Ny+Hy] = FT.(grid_info.RAC)
    hFW[1+Hx:Nx+Hx, 1+Hy:Ny+Hy, 1+Hz:Nz+Hz] = FT.(grid_info.hFacW)
    hFS[1+Hx:Nx+Hx, 1+Hy:Ny+Hy, 1+Hz:Nz+Hz] = FT.(grid_info.hFacS)
    hFC[1+Hx:Nx+Hx, 1+Hy:Ny+Hy, 1+Hz:Nz+Hz] = FT.(grid_info.hFacC)

    ##### fill halos
    @views @. dzF[1:Hz] = dzF[Hz+1]
    @views @. dzC[1:Hz] = dzC[Hz+1]
    @views @. dzF[Nz+Hz+1:Nz+2*Hz] = dzF[Nz+Hz]
    @views @. dzC[Nz+Hz+1:Nz+2*Hz] = dzC[Nz+Hz]

    @views @. zF[1:Hz] = zF[1+Hz:2*Hz] + dzF[1]*2 
    @views @. zC[1:Hz] = zC[1+Hz:2*Hz] + dzC[1]*2
    @views @. zF[Nz+Hz+2:Nz+2*Hz] = zF[Nz+2:Nz+Hz] - dzF[end]*2
    @views @. zC[Nz+Hz+1:Nz+2*Hz] = zC[Nz+1:Nz+Hz] - dzC[end]*2

    if TX == Periodic
        for g in (dxC, dyC, dxF, dyF, Az)
            @views @. g[1:Hx, :] = g[Nx+1:Nx+Hx, :]            # west
            @views @. g[Nx+Hx+1:Nx+2Hx, :] = g[1+Hx:2Hx, :]    # east
        end
    elseif TX == Bounded
        for g in (dxC, dyC, dxF, dyF, Az)
            @views @. g[1:Hx, :] = g[Hx+1:Hx+1, :]             # west
            @views @. g[Nx+Hx+1:Nx+2Hx, :] = g[Nx+Hx:Nx+Hx, :] # east
        end
    end

    for g in (dxC, dyC, dxF, dyF, Az)
        @views @. g[:, 1:Hy] = g[:, Hy+1:Hy+1]             # south
        @views @. g[:, Ny+Hy+1:Ny+2Hy] = g[:, Ny+Hy:Ny+Hy] # north
    end

    for g in (hFW, hFS, hFC)
        if TX == Periodic
            @views @. g[1:Hx, :, :] = g[Nx+1:Nx+Hx, :, :]            # west
            @views @. g[Nx+Hx+1:Nx+2Hx, :, :] = g[1+Hx:2Hx, :, :]    # east
        elseif TX == Bounded
            @views @. g[1:Hx, :, :] = g[Hx+1:Hx+1, :, :]             # west
            @views @. g[Nx+Hx+1:Nx+2Hx, :, :] = g[Nx+Hx:Nx+Hx, :, :] # east
        end
        @views @. g[:, 1:Hy, :] = g[:, Hy+1:Hy+1, :]             # south
        @views @. g[:, Ny+Hy+1:Ny+2Hy, :] = g[:, Ny+Hy:Ny+Hy, :] # north
        @views @. g[:, :, 1:Hz] = g[:, :, Hz+1:Hz+1]             # top
        @views @. g[:, :, Nz+Hz+1:Nz+2Hz] = g[:, :, Nz+Hz:Nz+Hz] # bottom
    end

    ##### calculate areas and volumes
    for i in 1:Nx+2*Hx
        for j in 1:Ny+2*Hy
            for k in 1:Nz+2*Hz
                Ax[i,j,k] = dzF[k] * dyF[i,j] * hFW[i,j,k]
                Ay[i,j,k] = dzF[k] * dxF[i,j] * hFS[i,j,k]
                Vol[i,j,k] = dzF[k] * Az[i,j] * hFC[i,j,k]
            end
        end
    end

    landmask = landmask_validation(landmask, Nx, Ny, Nz, Hx, Hy, Hz, FT, TX, TY)

    return LatLonGrid{FT, TX, TY, TZ, typeof(xF), typeof(dxC), typeof(Vol)}(
        xC, yC, zC, xF, yF, zF, Δx, Δy, dxC, dyC, dzC, dxF, dyF, dzF, Ax, Ay, Az, Vol, Nx, Ny, Nz, Hx, Hy, Hz, landmask)
end

function show(io::IO, g::LatLonGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    xL, xR = g.xF[g.Hx+1], g.xF[g.Hx+1+g.Nx]
    yL, yR = g.yF[g.Hy+1], g.yF[g.Hy+1+g.Ny]
    zL, zR = g.zF[g.Hz+1], g.zF[g.Hz+1+g.Nz]
    dzF_min = minimum(g.dzF)
    dzF_max = maximum(g.dzF)
    print(io, "LatLonGrid{$FT, $TX, $TY, $TZ}\n",
              "domain: x ∈ [$xL, $xR], y ∈ [$yL, $yR], z ∈ [$zL, $zR]\n",
              "topology (Tx, Ty, Tz):     ", (TX, TY, TZ), '\n',
              "resolution (Nx, Ny, Nz):   ", (g.Nx, g.Ny, g.Nz), '\n',
              "halo size (Hx, Hy, Hz):    ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): ", g.Δx, ", ", g.Δy, ", [min=", dzF_min, ", max=", dzF_max,"])")
end

function short_show(grid::LatLonGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    return "LatLonGrid{$FT, $TX, $TY, $TZ}(Nx=$(grid.Nx), Ny=$(grid.Ny), Nz=$(grid.Nz))"
end
