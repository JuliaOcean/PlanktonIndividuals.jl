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
