const R_Earth = 6371.0e3    # Mean radius of the Earth in meter (https://en.wikipedia.org/wiki/Earth)

struct RegularLatLonGrid{TX, TY, TZ, R} <: AbstractGrid{TX, TY, TZ}
    # corrdinates at cell centers, unit: degree
    xC::R
    yC::R
    # corrdinates at cell centers, unit: meter
    zC::R
    # corrdinates at cell faces, unit: degree
    xF::R
    yF::R
    # corrdinates at cell faces, unit: meter
    zF::R
    # grid spacing, unit: degree
    Δx::Float64
    Δy::Float64
    # grid spacing, unit: meter
    Δz::Float64
    # number of grid points
    Nx::Int
    Ny::Int
    Nz::Int
    # number of halo points
    Hx::Int
    Hy::Int
    Hz::Int
    # radius, unit: meter
    radius::Float64
end

"""
    RegularLatLonGrid(;size, lat, lon,
                       radius = R_Earth,
                       halo = (2, 2, 2))
Creats a `RegularLatLonGrid` struct with `size = (Nx, Ny, Nz)` grid points.

Keyword Arguments (Required)
============================
- `size` : A tuple prescribing the number of grid points. 
                `size` is a 3-tuple no matter for 3D, 2D, or 1D model.
- `lat` : A 2-tuple specifying the startind and ending points in latitudinal direction.
                Possible values are from -80 (80S) to 80 (80N).
- `lon` : A 2-tuple specifying the startind and ending points in longitudinal direction.
                Possible values are from -180 (180W) to 180 (180E).
- `z` : A 2-tuple specifying the startind and ending points in vertical direction.
                Vertical indexing starts from surface and use negative numbers for depth.

Keyword Arguments (Optional)
============================
- `radius` : Specify the radius of the Earth used in the model.
- `halo` : A tuple of integers that specifies the size of the halo region of cells
                surrounding the physical interior for each direction.
                `halo` is a 3-tuple no matter for 3D, 2D, or 1D model.
                At least 2 halo points are needed for DST3FL advection scheme.
"""
function RegularLatLonGrid(;size, lat, lon, z,
                            radius = R_Earth,
                            halo = (2, 2, 2))
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    lat₁, lat₂ = lat
    lon₁, lon₂ = lon
    z₁, z₂ = z

    @assert -180 <= lon₁ < lon₂ <= 180
    @assert -80 <= lat₁ < lat₂ <= 80
    @assert z₁ > z₂

    TX = lon₁ == -180 && lon₂ == 180 ? Periodic : Bounded
    TY = Bounded
    TZ = Bounded

    Δx = (lon₂ - lon₁) / Nx
    Δy = (lat₂ - lat₁) / Ny
    Δz = -(z₂ - z₁) / Nz

    xF = range(lon₁ - Hx * Δx, lon₁ + (Nx + Hx - 1) * Δx, length = Nx + 2 * Hx)
    yF = range(lat₁ - Hy * Δy, lat₁ + (Ny + Hy - 1) * Δy, length = Ny + 2 * Hy)
    zF = -range(z₁ - Hz * Δz, z₂ + (Nz + Hz - 1) * Δz, length = Nz + 2 * Hz)

    xC = range(lon₁ + (0.5 - Hx) * Δx, lon₁ + (Nx + Hx - 0.5) * Δx, length = Nx + 2 * Hx)
    yC = range(lat₁ + (0.5 - Hy) * Δy, lat₁ + (Ny + Hy - 0.5) * Δy, length = Ny + 2 * Hy)
    zC = -range(z₁ + (0.5 - Hz) * Δz, z₂ + (Nz + Hz - 0.5) * Δz, length = Nz + 2 * Hz)

    return RegularLatLonGrid{TX, TY, TZ, typeof(xF)}(
        xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Nx, Ny, Nz, Hx, Hy, Hz, radius)
end

import Base: show

function show(io::IO, g::RegularLatLonGrid{TX, TY, TZ}) where {TX, TY, TZ}
    xL, xR = g.xF[g.Hx+1], g.xF[g.Hx+1+g.Nx]
    yL, yR = g.yF[g.Hy+1], g.yF[g.Hy+1+g.Ny]
    zL, zR = g.zF[g.Hz+1], g.zF[g.Hz+1+g.Nz]
    print(io, "domain: x ∈ [$xL, $xR], y ∈ [$yL, $yR], z ∈ [$zL, $zR]\n",
              "topology (Tx, Ty, Tz):     ", (TX, TY, TZ), '\n',
              "resolution (Nx, Ny, Nz):   ", (g.Nx, g.Ny, g.Nz), '\n',
              "halo size (Hx, Hy, Hz):    ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): ", (g.Δx, g.Δy, g.Δz))
end