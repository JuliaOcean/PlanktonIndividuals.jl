
"""
    which_grid(phyt,g)
decide which grid an individual is in
return grid indices
"""
function which_grid(phyt, g)
    x = phyt[1]; y = phyt[2]; z = phyt[3];
    xind = findall(t -> t≤x, g.xF[2:end-1])[end]
    yind = findall(t -> t≤y, g.yF[2:end-1])[end]
    zind = findall(t -> t≤z, g.zF[2:end-2])[end]
    return xind, yind, zind
end

"""
    count_chl(phyts_a, grid)
Compute total chl concentration in each grid cell
"""
function count_chl(phyts_a, grid)
    cells = zeros(grid.Nx, grid.Ny, grid.Nz)
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        x,y,z = which_grid(phyt, grid)
        cells[x, y, z] = cells[x, y, z] + phyt[9]
    end
    cells .= cells ./ grid.V
    return cells
end

"""
    count_vertical_num(phyts_a)
Count individual numbers in each meter vertically, accumulate horizontally
"""
function count_vertical_num(phyts_a)
    VD = zeros(500)
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        z = -trunc(Int, phyt[3]) + 1
        VD[z] = VD[z] + 1.0
    end
    return VD
end

"""
    count_horizontal_num(phyts_a, grid)
Count individual numbers in each horizontal grid cell, accumulate in vertical direction
"""
function count_horizontal_num(phyts_a,grid)
    HD = zeros(grid.Nx,grid.Ny)
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        x,y,z = which_grid(phyt, gird)
        HD[x,y] = HD[x,y] + 1.0
    end
    return HD
end
