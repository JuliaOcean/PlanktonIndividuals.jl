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
@kernel function periodic_domain⁺_kernel!(op_array, g::Grids, ind)
    i = @index(Global, Linear)
    if op_array[i,1+ind] ≥ g.xF[g.Nx+g.Hx+1]
        op_array[i,1+ind] = op_array[i,1+ind] - g.Nx*g.Δx
        # op_array[i,1+ind] = g.xF[g.Nx+g.Hx+1]
    end
    if op_array[i,2+ind] ≥ g.yF[g.Ny+g.Hy+1]
        op_array[i,2+ind] = op_array[i,2+ind] - g.Ny*g.Δy
        # op_array[i,2+ind] = g.yF[g.Ny+g.Hy+1]
    end
    if op_array[i,3+ind] ≥ g.zF[g.Nz+g.Hz+1]
        op_array[i,3+ind] = g.zF[g.Nz+g.Hz+1]
    end
    # @inbounds op_array[i,3+ind] = min(g.zF[g.Hz+g.Nz+1], op_array[i,3+ind])
end
function periodic_domain⁺!(op_array, arch::Architecture, g::Grids, ind::Int64)
    kernel! = periodic_domain⁺_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, g, ind)
    wait(device(arch), event)
    return nothing
end
@kernel function periodic_domain⁻_kernel!(op_array, g::Grids, ind)
    i = @index(Global, Linear)
    if op_array[i,1+ind] < g.xF[g.Hx+1]
        # op_array[i,1+ind] = g.xF[g.Hx+1]
        op_array[i,1+ind] = op_array[i,1+ind] + g.Nx*g.Δx
    end
    if op_array[i,2+ind] < g.yF[g.Hy+1]
        # op_array[i,2+ind] = g.yF[g.Hy+1]
        op_array[i,2+ind] = op_array[i,2+ind] + g.Ny*g.Δy
    end
    if op_array[i,3+ind] < g.zF[g.Hz+1]
        op_array[i,3+ind] = g.zF[g.Hz+1]
    end
    # @inbounds op_array[i,1+ind] = op_array[i,1+ind] > g.xF[g.Hx+1] ? op_array[i,1+ind] : op_array[i,1+ind] + g.Nx*g.Δx
    # @inbounds op_array[i,2+ind] = op_array[i,2+ind] > g.yF[g.Hy+1] ? op_array[i,2+ind] : op_array[i,2+ind] + g.Ny*g.Δy
    # @inbounds op_array[i,3+ind] = max(g.zF[g.Hz+1], op_array[i,3+ind])
end
function periodic_domain⁻!(op_array, arch::Architecture, g::Grids, ind::Int64)
    kernel! = periodic_domain⁻_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, g, ind)
    wait(device(arch), event)
    return nothing
end
function in_domain!(op_array, arch::Architecture, g::Grids, ind)
    periodic_domain⁺!(op_array, arch::Architecture, g::Grids, ind)
    periodic_domain⁻!(op_array, arch::Architecture, g::Grids, ind)
end

in_domain!(op_array, arch::Architecture, g::Grids) = in_domain!(op_array, arch::Architecture, g::Grids, 0)


##### find indices (halo points excluded)
@inline find_xF_ind(x, g::Grids) = x ≥ g.xF[g.Hx+1] ? x ÷ g.Δx + 1 : 0.0

@inline find_yF_ind(y, g::Grids) = y ≥ g.yF[g.Hy+1] ? y ÷ g.Δy + 1 : 0.0

@inline find_zF_ind(z, g::Grids) = z ≥ g.zF[g.Hz+1] ? (g.Nz * g.Δz + z) ÷ g.Δz + 1 : 0.0

@inline find_xC_ind(x, g::Grids) = x ≥ g.xC[g.Hx+1] ? (x - g.Δx * 0.5) ÷ g.Δx + 1 : 0.0

@inline find_yC_ind(y, g::Grids) = y ≥ g.yC[g.Hy+1] ? (y - g.Δy * 0.5) ÷ g.Δy + 1 : 0.0

@inline find_zC_ind(z, g::Grids) = z ≥ g.zC[g.Hz+1] ? (g.Nz * g.Δz + z - g.Δz * 0.5) ÷ g.Δz + 1 : 0.0

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
    @inbounds op_array[i,12] = (op_array[i,1] - g.xF[x₀]) / g.Δx # dx
    @inbounds op_array[i,13] = (op_array[i,2] - g.yC[y₀]) / g.Δy # dy
    @inbounds op_array[i,14] = (op_array[i,3] - g.zC[z₀]) / g.Δz # dz
end
@kernel function find_yᵈ_kernel!(op_array, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds op_array[i,12] = (op_array[i,1] - g.xC[x₀]) / g.Δx # dx
    @inbounds op_array[i,13] = (op_array[i,2] - g.yF[y₀]) / g.Δy # dy
    @inbounds op_array[i,14] = (op_array[i,3] - g.zC[z₀]) / g.Δz # dz
end
@kernel function find_zᵈ_kernel!(op_array, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds op_array[i,12] = (op_array[i,1] - g.xC[x₀]) / g.Δx # dx
    @inbounds op_array[i,13] = (op_array[i,2] - g.yC[y₀]) / g.Δy # dy
    @inbounds op_array[i,14] = (op_array[i,3] - g.zF[z₀]) / g.Δz # dz
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

##### velocity interpolation for each individual
@inline trilinear_itpl(u000, u100, u010, u110, u001, u101, u011, u111, dx, dy, dz) =
    ((u000 * (1 - dx) + u100 * dx) * (1 - dy) + (u010 * (1-dx) + u110 * dx) * dy) * (1 - dz) +
    ((u001 * (1 - dx) + u101 * dx) * (1 - dy) + (u011 * (1-dx) + u111 * dx) * dy) * dz

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
@kernel function calc_intermediate_coord_kernel!(op_array, vel_array, ΔT, ind::Int64, dt::Float64)
    i = @index(Global, Linear)
    @inbounds op_array[i,15] = op_array[i,1] + dt * vel_array[i,1+ind*3] * ΔT
    @inbounds op_array[i,16] = op_array[i,2] + dt * vel_array[i,2+ind*3] * ΔT
    @inbounds op_array[i,17] = op_array[i,3] + dt * vel_array[i,3+ind*3] * ΔT
end
function calc_intermediate_coord!(op_array, vel_array, arch::Architecture, ΔT, ind::Int64, dt::Float64)
    kernel! = calc_intermediate_coord_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, vel_array, ΔT, ind, dt)
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
@kernel function calc_coord_kernel!(phytos, vel_array, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[i,1] = phytos[i,1] + vel_array[i,1] * ΔT
    @inbounds phytos[i,2] = phytos[i,2] + vel_array[i,2] * ΔT
    @inbounds phytos[i,3] = phytos[i,3] + vel_array[i,3] * ΔT
end
function calc_coord!(phytos, vel_array, arch::Architecture, ΔT)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(vel_array,1),))
    event = kernel!(phytos, vel_array, ΔT)
    wait(device(arch), event)
    return nothing
end
