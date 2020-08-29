"""
    grid_Ogrids(Ogrid, filepath)
Read grid information from Oceananigans
Return a grid 'struc'
'z' starts from the bottom
write grid information to a netCDF file
"""
function read_Ogrids(Ogrid)
    Nx = Ogrid.Nx
    Ny = Ogrid.Ny
    Nz = Ogrid.Nz
    Hx = Ogrid.Hx
    Hy = Ogrid.Hy
    Hz = Ogrid.Hz
    xC = collect(Ogrid.xC)
    yC = collect(Ogrid.yC)
    zC = collect(Ogrid.zC)
    xF = collect(Ogrid.xF)
    yF = collect(Ogrid.yF)
    zF = collect(Ogrid.zF)
    Δx = Ogrid.Δx
    Δy = Ogrid.Δy
    Δz = Ogrid.Δz
    Ax = Δy*Δz
    Ay = Δx*Δz
    Az = Δx*Δy
    V  = Δx*Δy*Δz
    g = grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Ax, Ay, Az, V, Nx, Ny, Nz, Hx, Hy, Hz)
    return g
end
function read_Ogrids(Ogrid, filepath)
    Nx = Ogrid.Nx
    Ny = Ogrid.Ny
    Nz = Ogrid.Nz
    Hx = Ogrid.Hx
    Hy = Ogrid.Hy
    Hz = Ogrid.Hz
    xC = collect(Ogrid.xC)
    yC = collect(Ogrid.yC)
    zC = collect(Ogrid.zC)
    xF = collect(Ogrid.xF)
    yF = collect(Ogrid.yF)
    zF = collect(Ogrid.zF)
    Δx = Ogrid.Δx
    Δy = Ogrid.Δy
    Δz = Ogrid.Δz
    Ax = Δy*Δz
    Ay = Δx*Δz
    Az = Δx*Δy
    V  = Δx*Δy*Δz

    ds = NCDataset(filepath*"grid.nc","c")
    ds.attrib["Nx"] = Nx
    ds.attrib["Ny"] = Ny
    ds.attrib["Nz"] = Nz
    ds.attrib["Hx"] = Hx
    ds.attrib["Hy"] = Hy
    ds.attrib["Hz"] = Hz
    ds.attrib["dx"] = Δx
    ds.attrib["dy"] = Δy
    ds.attrib["dz"] = Δz
    defVar(ds,"xC",xC,("xC",), attrib = Dict("units" => "m"))
    defVar(ds,"yC",yC,("yC",), attrib = Dict("units" => "m"))
    defVar(ds,"zC",zC,("zC",), attrib = Dict("units" => "m"))
    defVar(ds,"xF",xF,("xF",), attrib = Dict("units" => "m"))
    defVar(ds,"yF",yF,("yF",), attrib = Dict("units" => "m"))
    defVar(ds,"zF",zF,("zF",), attrib = Dict("units" => "m"))
    close(ds)

    g = grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Ax, Ay, Az, V, Nx, Ny, Nz, Hx, Hy, Hz)
    return g
end

function gen_Grid(;size, spacing ,halo = (1, 1, 1),)
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    Δx, Δy, Δz = spacing

    xF = range(-Hx * Δx, Nx * Δx, length = Nx + 2 * Hx)
    yF = range(-Hy * Δy, Ny * Δy, length = Ny + 2 * Hy)
    zF = range(-(Hz + Nz) * Δz, Δz, length = Nz + 2 * Hz + 1)

    xC = range((0.5 - Hx) * Δx, (Nx + 0.5) * Δx, length = Nx + 2 * Hx)
    yC = range((0.5 - Hy) * Δy, (Ny + 0.5) * Δy, length = Ny + 2 * Hy)
    zC = range((0.5 - Hz - Nz) * Δz, 0.5 * Δz, length = Nz + 2 * Hz)

    Ax = Δy*Δz
    Ay = Δx*Δz
    Az = Δx*Δy
    V  = Δx*Δy*Δz

    return grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Ax, Ay, Az, V, Nx, Ny, Nz, Hx, Hy, Hz)
end

import Base: show

function show(io::IO, g::grids)
    xL, xR = g.xF[2], g.xF[end]
    yL, yR = g.yF[2], g.yF[end]
    zL, zR = g.zF[2], g.zF[end-1]
    print(io, "domain: x ∈ [$xL, $xR], y ∈ [$yL, $yR], z ∈ [$zL, $zR]\n",
              "resolution (Nx, Ny, Nz):   ", (g.Nx, g.Ny, g.Nz), '\n',
              "halo size (Hx, Hy, Hz):    ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): ", (g.Δx, g.Δy, g.Δz))
end
