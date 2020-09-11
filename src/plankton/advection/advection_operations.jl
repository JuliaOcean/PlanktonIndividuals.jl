##### set up the operating array(cuarray) for plankton advection
#=
index of op_array
1  2  3  4   5   6   7   8   9   10  11  12  13  14  15  16  17  18
x  y  z  xF  yF  zF  u0  u1  v0  v1  w0  w1  xd  yd  zd  xt  yt  zt

19  20  21  22  23  24  25  26  27  28  29  30
u1  v1  w1  u2  v2  w2  u3  v3  w3  u4  v4  w4
=#
@kernel function update_op_array_kernel!(phytos, ope)
    i = @index(Global, Linear)
    ope[i,1] = phytos[i,1]
    ope[i,2] = phytos[i,2]
    ope[i,3] = phytos[i,3]
end
function update_op_array!(phytos, ope, arch::Architecture)
    kernel! = update_op_array_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(phytos,ope)
    wait(device(arch), event)
    return nothing
end

##### deal with particles moved out of the domain
@kernel function periodic_domain⁺_kernel!(ope, g::Grids, ind)
    i = @index(Global, Linear)
    if ope[i,1+ind] ≥ g.xF[g.Nx+g.Hx+1]
        ope[i,1+ind] = ope[i,1+ind] - g.Nx*g.Δx
        # ope[i,1+ind] = g.xF[g.Nx+g.Hx+1]
    end
    if ope[i,2+ind] ≥ g.yF[g.Ny+g.Hy+1]
        ope[i,2+ind] = ope[i,2+ind] - g.Ny*g.Δy
        # ope[i,2+ind] = g.yF[g.Ny+g.Hy+1]
    end
    if ope[i,3+ind] ≥ g.zF[g.Nz+g.Hz+1]
        ope[i,3+ind] = g.zF[g.Nz+g.Hz+1]
    end
end
function periodic_domain⁺!(ope, arch::Architecture, g::Grids, ind::Int64)
    kernel! = periodic_domain⁺_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, g, ind)
    wait(device(arch), event)
    return nothing
end
@kernel function periodic_domain⁻_kernel!(ope, g::Grids, ind)
    i = @index(Global, Linear)
    if ope[i,1+ind] < g.xF[g.Hx+1]
        # ope[i,1+ind] = g.xF[g.Hx+1]
        ope[i,1+ind] = ope[i,1+ind] + g.Nx*g.Δx
    end
    if ope[i,2+ind] < g.yF[g.Hy+1]
        # ope[i,2+ind] = g.yF[g.Hy+1]
        ope[i,2+ind] = ope[i,2+ind] + g.Ny*g.Δy
    end
    if ope[i,3+ind] < g.zF[g.Hz+1]
        ope[i,3+ind] = g.zF[g.Hz+1]
    end
end
function periodic_domain⁻!(ope, arch::Architecture, g::Grids, ind::Int64)
    kernel! = periodic_domain⁻_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, g, ind)
    wait(device(arch), event)
    return nothing
end
function in_domain!(ope, arch::Architecture, g::Grids, ind)
    periodic_domain⁺!(ope, arch::Architecture, g::Grids, ind)
    periodic_domain⁻!(ope, arch::Architecture, g::Grids, ind)
end

in_domain!(ope, arch::Architecture, g::Grids) = in_domain!(ope, arch::Architecture, g::Grids, 0)


##### find indices (halo points excluded)
@inline find_xF_ind(x, g::Grids) = x ≥ g.xF[g.Hx+1] ? x ÷ g.Δx + 1 : 0.0

@inline find_yF_ind(y, g::Grids) = y ≥ g.yF[g.Hy+1] ? y ÷ g.Δy + 1 : 0.0

@inline find_zF_ind(z, g::Grids) = z ≥ g.zF[g.Hz+1] ? (g.Nz * g.Δz + z) ÷ g.Δz + 1 : 0.0

##### find cell and face indices for each individual
@kernel function find_inds_kernel!(ope, g::Grids)
    i = @index(Global, Linear)
    @inbounds ope[i,4] = find_xF_ind(ope[i,1], g)  # xF index
    @inbounds ope[i,5] = find_yF_ind(ope[i,2], g)  # yF index
    @inbounds ope[i,6] = find_zF_ind(ope[i,3], g)  # zF index
end
function find_inds!(ope, arch::Architecture, g::Grids)
    kernel! = find_inds_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, g)
    wait(device(arch), event)
    return nothing
end
@kernel function update_inds_kernel!(ope, g::Grids)
    i = @index(Global, Linear)
    @inbounds ope[i,4] = find_xF_ind(ope[i,16], g)  # xF index
    @inbounds ope[i,5] = find_yF_ind(ope[i,17], g)  # yF index
    @inbounds ope[i,6] = find_zF_ind(ope[i,18], g)  # zF index
end
function update_inds!(ope, arch::Architecture, g::Grids)
    kernel! = update_inds_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, g)
    wait(device(arch), event)
    return nothing
end

##### find velocities around the individual to interpolate
@kernel function find_vel_kernel!(ope, ind_array::AbstractArray{Int64,2}, g::Grids, u, v, w)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds ope[i,7]  = u[x₀,   y₀,   z₀  ]
    @inbounds ope[i,8]  = u[x₀+1, y₀,   z₀  ]
    @inbounds ope[i,9]  = v[x₀,   y₀,   z₀  ]
    @inbounds ope[i,10] = v[x₀,   y₀+1, z₀  ]
    @inbounds ope[i,11] = w[x₀,   y₀,   z₀  ]
    @inbounds ope[i,12] = w[x₀,   y₀,   z₀+1]
end
function find_vel!(ope, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids, u, v, w)
    kernel! = find_vel_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, ind_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

@kernel function find_xᵈ_kernel!(ope, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds ope[i,13] = (ope[i,1] - g.xF[x₀]) / g.Δx
    @inbounds ope[i,14] = (ope[i,2] - g.yF[y₀]) / g.Δy
    @inbounds ope[i,15] = (ope[i,3] - g.zF[z₀]) / g.Δz
end
function find_xᵈ!(ope, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids)
    kernel! = find_xᵈ_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, ind_array, g)
    wait(device(arch), event)
    return nothing
end

##### velocity interpolation for each individual
@inline linear_itpl(u0, u1, xd) = u0 * (1.0 - xd) + u1 * xd

@kernel function vel_interpolation_kernel!(ope, g::Grids, ind::Int64)
    i = @index(Global, Linear)
    @inbounds ope[i, 19+ind*3] = linear_itpl(ope[i,7],  ope[i,8],  ope[i,13])
    @inbounds ope[i, 20+ind*3] = linear_itpl(ope[i,9],  ope[i,10], ope[i,14])
    @inbounds ope[i, 21+ind*3] = linear_itpl(ope[i,11], ope[i,12], ope[i,15])
end
function vel_interpolation!(ope, arch::Architecture, g::Grids, ind::Int64)
    kernel! = vel_interpolation_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, g, ind)
    wait(device(arch), event)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_1st_intermediate_coord_kernel!(ope, ΔT)
    i = @index(Global, Linear)
    @inbounds ope[i,16] = ope[i,1] + 0.5 * ope[i,19] * ΔT
    @inbounds ope[i,17] = ope[i,2] + 0.5 * ope[i,20] * ΔT
    @inbounds ope[i,18] = ope[i,3] + 0.5 * ope[i,21] * ΔT
end
@kernel function calc_2nd_intermediate_coord_kernel!(ope, ΔT)
    i = @index(Global, Linear)
    @inbounds ope[i,16] = ope[i,16] + 0.5 * ope[i,22] * ΔT
    @inbounds ope[i,17] = ope[i,17] + 0.5 * ope[i,23] * ΔT
    @inbounds ope[i,18] = ope[i,18] + 0.5 * ope[i,24] * ΔT
end
@kernel function calc_3rd_intermediate_coord_kernel!(ope, ΔT)
    i = @index(Global, Linear)
    @inbounds ope[i,16] = ope[i,16] + 1.0 * ope[i,25] * ΔT
    @inbounds ope[i,17] = ope[i,17] + 1.0 * ope[i,26] * ΔT
    @inbounds ope[i,18] = ope[i,18] + 1.0 * ope[i,27] * ΔT
end
function calc_1st_intermediate_coord!(ope, arch::Architecture, ΔT)
    kernel! = calc_1st_intermediate_coord_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, ΔT)
    wait(device(arch), event)
    return nothing
end
function calc_2nd_intermediate_coord!(ope, arch::Architecture, ΔT)
    kernel! = calc_2nd_intermediate_coord_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, ΔT)
    wait(device(arch), event)
    return nothing
end
function calc_3rd_intermediate_coord!(ope, arch::Architecture, ΔT)
    kernel! = calc_3rd_intermediate_coord_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate final velocities by RK4
@kernel function calc_vel_rk4_kernel!(ope)
    i = @index(Global, Linear)
    @inbounds ope[i,19] = (ope[i,19] + 2*ope[i,22] + 2*ope[i,25] + ope[i,28]) /6
    @inbounds ope[i,20] = (ope[i,20] + 2*ope[i,23] + 2*ope[i,26] + ope[i,29]) /6
    @inbounds ope[i,21] = (ope[i,21] + 2*ope[i,24] + 2*ope[i,27] + ope[i,30]) /6
end
function calc_vel_rk4!(ope, arch::Architecture)
    kernel! = calc_vel_rk4_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope)
    wait(device(arch), event)
    return nothing
end

##### calculate coordinates of each individual
@kernel function calc_coord_kernel!(phytos, ope, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[i,1] = phytos[i,1] + ope[i,19] * ΔT
    @inbounds phytos[i,2] = phytos[i,2] + ope[i,20] * ΔT
    @inbounds phytos[i,3] = phytos[i,3] + ope[i,21] * ΔT
end
function calc_coord!(phytos, ope, arch::Architecture, ΔT)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(phytos,1),))
    event = kernel!(phytos, ope, ΔT)
    wait(device(arch), event)
    return nothing
end
