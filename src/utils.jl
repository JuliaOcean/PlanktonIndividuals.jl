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

function vel_copy!(vel::NamedTuple, O_vels, arch::Architecture)
    u = O_vels.u.data.parent |> array_type(arch)
    v = O_vels.v.data.parent |> array_type(arch)
    w = O_vels.w.data.parent |> array_type(arch)
    vel.u.data .= u
    vel.v.data .= v
    vel.w.data .= w[:,:,1:end-1]
end

function vel_copy!(vel::NamedTuple, u, v, w, arch::Architecture, g::Grids)
    uvel = u |> array_type(arch)
    vvel = v |> array_type(arch)
    wvel = w |> array_type(arch)

    vel.u.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= uvel
    vel.v.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= vvel
    vel.w.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= wvel

    fill_halo_vel!(vel, g)
end
