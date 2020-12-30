abstract type AbstractGrid{FT, A} end
struct Grids{FT, A} <: AbstractGrid{FT, A}
    xC::A
    yC::A
    zC::A
    xF::A
    yF::A
    zF::A
    Δx::FT # unit: meter
    Δy::FT # unit: meter
    Δz::FT # unit: meter
    Ax::FT # unit: m²
    Ay::FT # unit: m²
    Az::FT # unit: m²
    V ::FT # unit: m³
    Nx::Int
    Ny::Int
    Nz::Int
    Hx::Int
    Hy::Int
    Hz::Int
end

function gen_Grid(;size, spacing, halo = (1, 1, 1),)
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    Δx, Δy, Δz = spacing

    xF = range(-Hx * Δx, (Nx + Hx - 1) * Δx, length = Nx + 2 * Hx)
    yF = range(-Hy * Δy, (Ny + Hy - 1) * Δy, length = Ny + 2 * Hy)
    zF = range(-(Hz + Nz) * Δz, (Hz -1) * Δz, length = Nz + 2 * Hz)

    xC = range((0.5 - Hx) * Δx, (Nx + Hx - 0.5) * Δx, length = Nx + 2 * Hx)
    yC = range((0.5 - Hy) * Δy, (Ny + Hy - 0.5) * Δy, length = Ny + 2 * Hy)
    zC = range((0.5 - Hz - Nz) * Δz, (Hz - 0.5) * Δz, length = Nz + 2 * Hz)

    Ax = Δy*Δz
    Ay = Δx*Δz
    Az = Δx*Δy
    V  = Δx*Δy*Δz

    return Grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Ax, Ay, Az, V, Nx, Ny, Nz, Hx, Hy, Hz)
end

# """
#     grid_Ogrids(Ogrid, filepath)
# Read grid information from Oceananigans
# Return a grid 'struc'
# 'z' starts from the bottom
# write grid information to a netCDF file
# """
# function read_Ogrids(Ogrid)
#     Nx = Ogrid.Nx
#     Ny = Ogrid.Ny
#     Nz = Ogrid.Nz
#     Hx = Ogrid.Hx
#     Hy = Ogrid.Hy
#     Hz = Ogrid.Hz
#     xC = collect(Ogrid.xC)
#     yC = collect(Ogrid.yC)
#     zC = collect(Ogrid.zC)
#     xF = collect(Ogrid.xF)
#     yF = collect(Ogrid.yF)
#     zF = collect(Ogrid.zF)
#     Δx = Ogrid.Δx
#     Δy = Ogrid.Δy
#     Δz = Ogrid.Δz
#     Ax = Δy*Δz
#     Ay = Δx*Δz
#     Az = Δx*Δy
#     V  = Δx*Δy*Δz
#     g = Grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Ax, Ay, Az, V, Nx, Ny, Nz, Hx, Hy, Hz)
#     return g
# end
# function read_Ogrids(Ogrid, filepath)
#     Nx = Ogrid.Nx
#     Ny = Ogrid.Ny
#     Nz = Ogrid.Nz
#     Hx = Ogrid.Hx
#     Hy = Ogrid.Hy
#     Hz = Ogrid.Hz
#     xC = collect(Ogrid.xC)
#     yC = collect(Ogrid.yC)
#     zC = collect(Ogrid.zC)
#     xF = collect(Ogrid.xF)
#     yF = collect(Ogrid.yF)
#     zF = collect(Ogrid.zF)
#     Δx = Ogrid.Δx
#     Δy = Ogrid.Δy
#     Δz = Ogrid.Δz
#     Ax = Δy*Δz
#     Ay = Δx*Δz
#     Az = Δx*Δy
#     V  = Δx*Δy*Δz

#     ds = NCDataset(filepath*"grid.nc","c")
#     ds.attrib["Nx"] = Nx
#     ds.attrib["Ny"] = Ny
#     ds.attrib["Nz"] = Nz
#     ds.attrib["Hx"] = Hx
#     ds.attrib["Hy"] = Hy
#     ds.attrib["Hz"] = Hz
#     ds.attrib["dx"] = Δx
#     ds.attrib["dy"] = Δy
#     ds.attrib["dz"] = Δz
#     defVar(ds,"xC",xC,("xC",), attrib = Dict("units" => "m"))
#     defVar(ds,"yC",yC,("yC",), attrib = Dict("units" => "m"))
#     defVar(ds,"zC",zC,("zC",), attrib = Dict("units" => "m"))
#     defVar(ds,"xF",xF,("xF",), attrib = Dict("units" => "m"))
#     defVar(ds,"yF",yF,("yF",), attrib = Dict("units" => "m"))
#     defVar(ds,"zF",zF,("zF",), attrib = Dict("units" => "m"))
#     close(ds)

#     g = Grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Ax, Ay, Az, V, Nx, Ny, Nz, Hx, Hy, Hz)
#     return g
# end

import Base: show

function show(io::IO, g::Grids)
    xL, xR = g.xF[g.Hx+1], g.xF[g.Hx+1+g.Nx]
    yL, yR = g.yF[g.Hy+1], g.yF[g.Hy+1+g.Ny]
    zL, zR = g.zF[g.Hz+1], g.zF[g.Hz+1+g.Nz]
    print(io, "domain: x ∈ [$xL, $xR], y ∈ [$yL, $yR], z ∈ [$zL, $zR]\n",
              "resolution (Nx, Ny, Nz):   ", (g.Nx, g.Ny, g.Nz), '\n',
              "halo size (Hx, Hy, Hz):    ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): ", (g.Δx, g.Δy, g.Δz))
end
