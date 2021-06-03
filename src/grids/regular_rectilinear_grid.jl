abstract type AbstractGrid{TX, TY, YZ} end
abstract type AbstractTopology end
struct Periodic <: AbstractTopology end
struct Bounded <: AbstractTopology end

struct RegularRectilinearGrid{TX, TY, TZ, R} <: AbstractGrid{TX, TY, TZ}
    # corrdinates at cell centers, unit: meter
    xC::R
    yC::R
    zC::R
    # corrdinates at cell faces, unit: meter
    xF::R
    yF::R
    zF::R
    # grid spacing, unit: meter
    Δx::Float64
    Δy::Float64
    Δz::Float64
    # number of grid points
    Nx::Int
    Ny::Int
    Nz::Int
    # number of halo points
    Hx::Int
    Hy::Int
    Hz::Int
end

"""
    RegularRectilinearGrid(;size, spacing,
                            topology = (Periodic, Periodic, Bounded),
                            halo = (2, 2, 2))
Creats a `RegularRectilinearGrid` struct with `size = (Nx, Ny, Nz)` grid points.

Keyword Arguments (Required)
============================
- `size` : A tuple prescribing the number of grid points. 
                `size` is a 3-tuple no matter for 3D, 2D, or 1D model.
- `spacing` : A tuple prescribing the length of each grid point in x, y, and z directions.
                `spacing` is a 3-tuple no matter for 3D, 2D, or 1D model.

Keyword Arguments (Optional)
============================
- `topology` : A 3-tuple specifying the topology of the domain.
                The topology can be either Periodic or Bounded in each direction.
- `halo` : A tuple of integers that specifies the size of the halo region of cells
                surrounding the physical interior for each direction.
                `halo` is a 3-tuple no matter for 3D, 2D, or 1D model.
                At least 2 halo points are needed for DST3FL advection scheme.
"""
function RegularRectilinearGrid(;size, spacing,
                                topology = (Periodic, Periodic, Bounded),
                                halo = (2, 2, 2))
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    Δx, Δy, Δz = spacing
    TX, TY, TZ = validate_topology(topology)

    xF = range(-Hx * Δx, (Nx + Hx - 1) * Δx, length = Nx + 2 * Hx)
    yF = range(-Hy * Δy, (Ny + Hy - 1) * Δy, length = Ny + 2 * Hy)
    zF = -range(-Hz * Δz, (Nz + Hz - 1) * Δz, length = Nz + 2 * Hz)

    xC = range((0.5 - Hx) * Δx, (Nx + Hx - 0.5) * Δx, length = Nx + 2 * Hx)
    yC = range((0.5 - Hy) * Δy, (Ny + Hy - 0.5) * Δy, length = Ny + 2 * Hy)
    zC = -range((0.5 - Hz) * Δz, (Nz + Hz - 0.5) * Δz, length = Nz + 2 * Hz)

    return RegularRectilinearGrid{TX, TY, TZ, typeof(xF)}(
        xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Nx, Ny, Nz, Hx, Hy, Hz)
end

import Base: show

function show(io::IO, g::RegularRectilinearGrid{TX, TY, TZ}) where {TX, TY, TZ}
    xL, xR = g.xF[g.Hx+1], g.xF[g.Hx+1+g.Nx]
    yL, yR = g.yF[g.Hy+1], g.yF[g.Hy+1+g.Ny]
    zL, zR = g.zF[g.Hz+1], g.zF[g.Hz+1+g.Nz]
    print(io, "domain: x ∈ [$xL, $xR], y ∈ [$yL, $yR], z ∈ [$zL, $zR]\n",
              "topology (Tx, Ty, Tz):     ", (TX, TY, TZ), '\n',
              "resolution (Nx, Ny, Nz):   ", (g.Nx, g.Ny, g.Nz), '\n',
              "halo size (Hx, Hy, Hz):    ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): ", (g.Δx, g.Δy, g.Δz))
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