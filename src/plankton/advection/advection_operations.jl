##### set up the operating array(cuarray) for plankton advection
#=
index of op_array
1  2  3  4      5     6     7     8     9     10    11    12  13  14  15  16  17
x  y  z  v000   v100  v010  v110  v001  v101  v011  v111  dx  dy  dz  xt  yt  zt

index of ind_array
1  2  3  4  5  6  7  8  9
xF yC zC xC yF zC xC yC zF

index of vel_array
1  2  3  4  5  6  7  8  9  10 11 12
u1 v1 w1 u2 v2 w2 u3 v3 w3 u4 v4 w4
=#
function adv_op_array_setup(phytos, arch::Architecture)
    total_num = size(phytos, 1)
    op_array = zeros(total_num, 17) |> array_type(arch)

    op_array[:, 1:3] .= phytos[:, 1:3]
    return op_array
end
function ind_array_setup(phytos, arch::Architecture)
    total_num = size(phytos, 1)
    ind_array = zeros(total_num, 9) |> array_type(arch)
    return ind_array
end
function vel_array_setup(phytos, arch::Architecture)
    total_num = size(phytos, 1)
    vel_array = zeros(total_num, 12) |> array_type(arch)
    return vel_array
end

##### deal with particles moved out of the domain
function periodic_domain_x(x, g::Grids)
    x = x ≤ g.xF[g.Hx+1]      ? x + g.xF[g.Nx+g.Hx+1] : x
    x = x ≥ g.xF[g.Nx+g.Hx+1] ? x - g.xF[g.Nx+g.Hx+1] : x
    return x
end
function periodic_domain_y(y, g::Grids)
    y = y ≤ g.yF[g.Hy+1]      ? y + g.yF[g.Ny+g.Hy+1] : y
    y = y ≥ g.yF[g.Ny+g.Hy+1] ? y - g.yF[g.Ny+g.Hy+1] : y
    return y
end
function bounded_domain_z(z, g::Grids)
    z = z ≤ g.zF[g.Hz+1]      ? g.zF[g.Hz+1] : z
    z = z ≥ g.zF[g.Nz+g.Hz+1] ? g.zF[g.Nz+g.Hz+1] : z
    return z
end

##### find indices (halo points excluded)
@inline find_xF_ind(x, g::Grids) = x ≥ 0.0 ? x ÷ g.Δx + 1 : x ÷ g.Δx

@inline find_yF_ind(y, g::Grids) = y ≥ 0.0 ? y ÷ g.Δy + 1 : y ÷ g.Δy

@inline find_zF_ind(z, g::Grids) = (g.Nz * g.Δz + z) ≥ 0.0 ?
    (g.Nz * g.Δz + z) ÷ g.Δz + 1 : (g.Nz * g.Δz + z) ÷ g.Δz

@inline find_xC_ind(x, g::Grids) = (x - g.Δx * 0.5) ≥ 0.0 ?
    (x - g.Δx * 0.5) ÷ g.Δx + 1 : (x - g.Δx * 0.5) ÷ g.Δx

@inline find_yC_ind(y, g::Grids) = (y - g.Δy * 0.5) ≥ 0.0 ?
    (y - g.Δy * 0.5) ÷ g.Δy + 1 : (y - g.Δy * 0.5) ÷ g.Δy

@inline find_zC_ind(z, g::Grids) = (g.Nz * g.Δz + z - g.Δz * 0.5) ≥ 0.0 ?
    (g.Nz * g.Δz + z - g.Δz * 0.5) ÷ g.Δz + 1 : (g.Nz * g.Δz + z - g.Δz * 0.5) ÷ g.Δz

##### find cell and face indices for each individual
@kernel function find_inds_kernel!(ind_array, op_array, g::Grids)
    i = @index(Global, Linear)
    @inbounds ind_array[i,1] = find_xF_ind(op_array[i,1], g) # xF index
    @inbounds ind_array[i,2] = find_yC_ind(op_array[i,2], g) # yC index
    @inbounds ind_array[i,3] = find_zC_ind(op_array[i,3], g) # zC index
    @inbounds ind_array[i,4] = find_xC_ind(op_array[i,1], g) # xC index
    @inbounds ind_array[i,5] = find_yF_ind(op_array[i,2], g) # yF index
    @inbounds ind_array[i,6] = find_zC_ind(op_array[i,3], g) # zC index
    @inbounds ind_array[i,7] = find_xC_ind(op_array[i,1], g) # xC index
    @inbounds ind_array[i,8] = find_yC_ind(op_array[i,2], g) # yC index
    @inbounds ind_array[i,9] = find_zF_ind(op_array[i,3], g) # zF index
end
function find_inds!(ind_array, op_array, arch::Architecture, g::Grids)
    kernel! = find_inds_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(ind_array, op_array, g)
    wait(device(arch), event)
    return nothing
end
@kernel function update_inds_kernel!(ind_array, op_array, g::Grids)
    i = @index(Global, Linear)
    @inbounds ind_array[i,1] = find_xF_ind(op_array[i,15], g)  # xF index
    @inbounds ind_array[i,2] = find_yC_ind(op_array[i,16], g)  # yC index
    @inbounds ind_array[i,3] = find_zC_ind(op_array[i,17], g)  # zC index
    @inbounds ind_array[i,4] = find_xC_ind(op_array[i,15], g)  # xC index
    @inbounds ind_array[i,5] = find_yF_ind(op_array[i,16], g)  # yF index
    @inbounds ind_array[i,6] = find_zC_ind(op_array[i,17], g)  # zC index
    @inbounds ind_array[i,7] = find_xC_ind(op_array[i,15], g)  # xC index
    @inbounds ind_array[i,8] = find_yC_ind(op_array[i,16], g)  # yC index
    @inbounds ind_array[i,9] = find_zF_ind(op_array[i,17], g)  # zF index
end
function update_inds!(ind_array, op_array, arch::Architecture, g::Grids)
    kernel! = update_inds_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(ind_array, op_array, g)
    wait(device(arch), event)
    return nothing
end

##### find velocities around the individual to interpolate
@kernel function find_vel_kernel!(op_array, ind_array::AbstractArray{Int64,2}, g::Grids, u)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz                     # xyz
    @inbounds op_array[i,4]  = u[x₀,   y₀,   z₀  ] # 000
    @inbounds op_array[i,5]  = u[x₀+1, y₀,   z₀  ] # 100
    @inbounds op_array[i,6]  = u[x₀,   y₀+1, z₀  ] # 010
    @inbounds op_array[i,7]  = u[x₀+1, y₀+1, z₀  ] # 110
    @inbounds op_array[i,8]  = u[x₀,   y₀,   z₀+1] # 001
    @inbounds op_array[i,9]  = u[x₀+1, y₀,   z₀+1] # 101
    @inbounds op_array[i,10] = u[x₀,   y₀+1, z₀+1] # 011
    @inbounds op_array[i,11] = u[x₀+1, y₀+1, z₀+1] # 111
end
function find_vel!(op_array, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids, u)
    kernel! = find_vel_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, ind_array, g, u)
    wait(device(arch), event)
    return nothing
end

@kernel function find_xᵈ_kernel!(op_array, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds op_array[i,12] = op_array[i,1] - g.xF[x₀] # dx
    @inbounds op_array[i,13] = op_array[i,2] - g.yC[y₀] # dy
    @inbounds op_array[i,14] = op_array[i,3] - g.zC[z₀] # dz
end
@kernel function find_yᵈ_kernel!(op_array, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds op_array[i,12] = op_array[i,1] - g.xC[x₀] # dx
    @inbounds op_array[i,13] = op_array[i,2] - g.yF[y₀] # dy
    @inbounds op_array[i,14] = op_array[i,3] - g.zC[z₀] # dz
end
@kernel function find_zᵈ_kernel!(op_array, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds op_array[i,12] = op_array[i,1] - g.xC[x₀] # dx
    @inbounds op_array[i,13] = op_array[i,2] - g.yC[y₀] # dy
    @inbounds op_array[i,14] = op_array[i,3] - g.zF[z₀] # dz
end
function find_xᵈ!(op_array, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids)
    kernel! = find_xᵈ_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, ind_array, g)
    wait(device(arch), event)
    return nothing
end
function find_yᵈ!(op_array, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids)
    kernel! = find_yᵈ_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, ind_array, g)
    wait(device(arch), event)
    return nothing
end
function find_zᵈ!(op_array, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids)
    kernel! = find_zᵈ_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, ind_array, g)
    wait(device(arch), event)
    return nothing
end

@inline trilinear_itpl(u000, u100, u010, u110, u001, u101, u011, u111, dx, dy, dz) =
    ((u000 * (1 - dx) + u100 * dx) * (1 - dy) + (u010 * (1-dx) + u110 * dx) * dy) * (1 - dz) +
    ((u001 * (1 - dx) + u101 * dx) * (1 - dy) + (u011 * (1-dx) + u111 * dx) * dy) * dz


##### velocity interpolation for each individual
@kernel function vel_interpolation_kernel!(vel_array, op_array, g::Grids, ind::Int64)
    i = @index(Global, Linear)
    @inbounds vel_array[i, ind] = trilinear_itpl(op_array[i,4], op_array[i,5], op_array[i,6], op_array[i,7],
                                                 op_array[i,8], op_array[i,9], op_array[i,10], op_array[i,11],
                                                 op_array[i,12], op_array[i,13], op_array[i,14])
end
function vel_interpolation!(vel_array, op_array, arch::Architecture, g::Grids, ind::Int64)
    kernel! = vel_interpolation_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(vel_array, op_array, g, ind)
    wait(device(arch), event)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_intermediate_coord_kernel!(op_array, vel_array, g::Grids, ΔT, ind::Int64)
    i = @index(Global, Linear)
    @inbounds op_array[i,15] = periodic_domain_x(op_array[i,1] + 0.5 * vel_array[i,1+ind*3] * ΔT, g)
    @inbounds op_array[i,16] = periodic_domain_y(op_array[i,2] + 0.5 * vel_array[i,2+ind*3] * ΔT, g)
    @inbounds op_array[i,17] =  bounded_domain_z(op_array[i,3] + 0.5 * vel_array[i,3+ind*3] * ΔT, g)
end
function calc_intermediate_coord!(op_array, vel_array, arch::Architecture, g::Grids, ΔT, ind::Int64)
    kernel! = calc_intermediate_coord_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, vel_array, g, ΔT, ind)
    wait(device(arch), event)
    return nothing
end

##### calculate final velocities by RK4
@kernel function calc_vel_rk4_kernel!(vel_array)
    i = @index(Global, Linear)
    @inbounds vel_array[i,1] = (vel_array[i,1] + 2*vel_array[i,4] + 2*vel_array[i,7] + vel_array[i,10]) /6
    @inbounds vel_array[i,2] = (vel_array[i,2] + 2*vel_array[i,5] + 2*vel_array[i,8] + vel_array[i,11]) /6
    @inbounds vel_array[i,3] = (vel_array[i,3] + 2*vel_array[i,6] + 2*vel_array[i,9] + vel_array[i,12]) /6
end
function calc_vel_rk4!(vel_array, arch::Architecture)
    kernel! = calc_vel_rk4_kernel!(device(arch), 256, (size(vel_array,1),))
    event = kernel!(vel_array)
    wait(device(arch), event)
    return nothing
end

##### calculate coordinates of each individual
@kernel function calc_coord_kernel!(phytos, vel_array, g::Grids, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[i,1] = periodic_domain_x(phytos[i,1] + vel_array[i,1] * ΔT, g)
    @inbounds phytos[i,2] = periodic_domain_y(phytos[i,2] + vel_array[i,2] * ΔT, g)
    @inbounds phytos[i,3] =  bounded_domain_z(phytos[i,3] + vel_array[i,3] * ΔT, g)
end
function calc_coord!(phytos, vel_array, arch::Architecture, g::Grids, ΔT)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(vel_array,1),))
    event = kernel!(phytos, vel_array, g, ΔT)
    wait(device(arch), event)
    return nothing
end
