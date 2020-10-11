##### copy external velocities into the model
function vel_copy!(vel::NamedTuple, u, v, w, arch::Architecture, g::Grids)
    vel.u.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= u
    vel.v.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= v
    vel.w.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= w

    fill_halo_vel!(vel, g)
end

"""
    sub_nut_tendency!(a, b, c)
subtract one tendency from total tendencies
"""
function sub_nut_tendency!(a, b, c)
    for i in 1:length(a)
        a[i].data .= b[i].data .- c[i].data
    end
end

function zero_fields!(a)
    for i in 1:length(a)
        a[i].data .= 0.0
    end
end
# """
#     add_nut_tendency!(a, b)
# add one tendency to another tendency
# """
# function add_nut_tendency!(a, b)
#     for i in 1:length(a)
#         a[i].data .= a[i].data .+ b[i].data
#     end
# end
