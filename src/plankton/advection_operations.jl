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

##### set up the operating array(cuarray) for plankton advection
function adv_op_array_setup(phytos, arch::Architecture)
    total_num = size(phytos, 2)
    op_array = zeros(24, total_num) |>array_type(arch)

    op_array[1:3, :] .= phytos[1:3, :]
    return op_array
end

##### find cell and face indices for each individual
@kernel function find_inds_kernel!(op_array, g::Grids)
    i = @index(Global, Linear)
    @inbounds op_array[4,i] = find_xF_ind(op_array[1,i], g)  # xF index
    @inbounds op_array[5,i] = find_xC_ind(op_array[1,i], g)  # xC index
    @inbounds op_array[6,i] = find_yF_ind(op_array[2,i], g)  # yF index
    @inbounds op_array[7,i] = find_yC_ind(op_array[2,i], g)  # yC index
    @inbounds op_array[8,i] = find_zF_ind(op_array[3,i], g)  # zF index
    @inbounds op_array[9,i] = find_zC_ind(op_array[3,i], g)  # zC index
end
function find_inds!(op_array, arch::Architecture, g::Grids)
    kernel! = find_inds_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g)
    wait(device(arch), event)
    return nothing
end
@kernel function update_inds_kernel!(op_array, g::Grids)
    i = @index(Global, Linear)
    @inbounds op_array[4,i] = find_xF_ind(op_array[10,i], g)  # xF index
    @inbounds op_array[5,i] = find_xC_ind(op_array[10,i], g)  # xC index
    @inbounds op_array[6,i] = find_yF_ind(op_array[11,i], g)  # yF index
    @inbounds op_array[7,i] = find_yC_ind(op_array[11,i], g)  # yC index
    @inbounds op_array[8,i] = find_zF_ind(op_array[12,i], g)  # zF index
    @inbounds op_array[9,i] = find_zC_ind(op_array[12,i], g)  # zC index
end
function update_inds!(op_array, arch::Architecture, g::Grids)
    kernel! = update_inds_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g)
    wait(device(arch), event)
    return nothing
end

##### velocity interpolation for each individual
@kernel function vel₁_interpolation_kernel!(op_array, g::Grids, u, v, w)
    i = @index(Global, Linear)
    @inbounds op_array[13,i] = trilinear_itlp(op_array[1:3,i], op_array[4,i], op_array[7,i], op_array[9,i], u)
    @inbounds op_array[14,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[6,i], op_array[9,i], v)
    @inbounds op_array[15,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[7,i], op_array[8,i], w)
end
function vel₁_interpolation!(op_array, arch::Architecture, g::Grids, u, v, w)
    kernel! = vel₁_interpolation_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

@kernel function vel₂_interpolation_kernel!(op_array, g::Grids, u, v, w)
    i = @index(Global, Linear)
    @inbounds op_array[16,i] = trilinear_itlp(op_array[1:3,i], op_array[4,i], op_array[7,i], op_array[9,i], u)
    @inbounds op_array[17,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[6,i], op_array[9,i], v)
    @inbounds op_array[18,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[7,i], op_array[8,i], w)
end
function vel₂_interpolation!(op_array, arch::Architecture, g::Grids, u, v, w)
    kernel! = vel₂_interpolation_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

@kernel function vel₃_interpolation_kernel!(op_array, g::Grids, u, v, w)
    i = @index(Global, Linear)
    @inbounds op_array[19,i] = trilinear_itlp(op_array[1:3,i], op_array[4,i], op_array[7,i], op_array[9,i], u)
    @inbounds op_array[20,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[6,i], op_array[9,i], v)
    @inbounds op_array[21,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[7,i], op_array[8,i], w)
end
function vel₃_interpolation!(op_array, arch::Architecture, g::Grids, u, v, w)
    kernel! = vel₃_interpolation_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

@kernel function vel₄_interpolation_kernel!(op_array, g::Grids, u, v, w)
    i = @index(Global, Linear)
    @inbounds op_array[22,i] = trilinear_itlp(op_array[1:3,i], op_array[4,i], op_array[7,i], op_array[9,i], u)
    @inbounds op_array[23,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[6,i], op_array[9,i], v)
    @inbounds op_array[24,i] = trilinear_itlp(op_array[1:3,i], op_array[5,i], op_array[7,i], op_array[8,i], w)
end
function vel₄_interpolation!(op_array, arch::Architecture, g::Grids, u, v, w)
    kernel! = vel₄_interpolation_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_intermediate1_coord_kernel!(op_array, g::Grids, ΔT)
    i = @index(Global, Linear)
    @inbounds op_array[10,i] = periodic_domain_x(op_array[1,i] + 0.5 * op_array[13,i] * ΔT, g)
    @inbounds op_array[11,i] = periodic_domain_y(op_array[2,i] + 0.5 * op_array[14,i] * ΔT, g)
    @inbounds op_array[12,i] =  bounded_domain_z(op_array[3,i] + 0.5 * op_array[15,i] * ΔT, g)
end
function calc_intermediate1_coord!(op_array, arch::Architecture, g::Grids, ΔT)
    kernel! = calc_intermediate1_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, ΔT)
    wait(device(arch), event)
    return nothing
end

@kernel function calc_intermediate2_coord_kernel!(op_array, g::Grids, ΔT)
    i = @index(Global, Linear)
    @inbounds op_array[10,i] = periodic_domain_x(op_array[1,i] + 0.5 * op_array[16,i] * ΔT, g)
    @inbounds op_array[11,i] = periodic_domain_y(op_array[2,i] + 0.5 * op_array[17,i] * ΔT, g)
    @inbounds op_array[12,i] =  bounded_domain_z(op_array[3,i] + 0.5 * op_array[18,i] * ΔT, g)
end
function calc_intermediate2_coord!(op_array, arch::Architecture, g::Grids, ΔT)
    kernel! = calc_intermediate2_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, ΔT)
    wait(device(arch), event)
    return nothing
end

@kernel function calc_intermediate3_coord_kernel!(op_array, g::Grids, ΔT)
    i = @index(Global, Linear)
    @inbounds op_array[10,i] = periodic_domain_x(op_array[1,i] + 0.5 * op_array[19,i] * ΔT, g)
    @inbounds op_array[11,i] = periodic_domain_y(op_array[2,i] + 0.5 * op_array[20,i] * ΔT, g)
    @inbounds op_array[12,i] =  bounded_domain_z(op_array[3,i] + 0.5 * op_array[21,i] * ΔT, g)
end
function calc_intermediate3_coord!(op_array, arch::Architecture, g::Grids, ΔT)
    kernel! = calc_intermediate3_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate final velocities by RK4
@kernel function calc_vel_rk4_kernel!(op_array)
    i = @index(Global, Linear)
    @inbounds op_array[13,i] = (op_array[13,i] + 2*op_array[16,i] + 2*op_array[19,i] + op_array[22,i]) /6
    @inbounds op_array[14,i] = (op_array[14,i] + 2*op_array[17,i] + 2*op_array[20,i] + op_array[23,i]) /6
    @inbounds op_array[15,i] = (op_array[15,i] + 2*op_array[18,i] + 2*op_array[21,i] + op_array[24,i]) /6
end
function calc_vel_rk4!(op_array, arch::Architecture)
    kernel! = calc_vel_rk4_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array)
    wait(device(arch), event)
    return nothing
end

##### calculate coordinates of each individual
@kernel function calc_coord_kernel!(op_array, phytos, g::Grids, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[1,i] = periodic_domain_x(op_array[1,i] + op_array[13,i] * ΔT, g)
    @inbounds phytos[2,i] = periodic_domain_y(op_array[2,i] + op_array[14,i] * ΔT, g)
    @inbounds phytos[3,i] =  bounded_domain_z(op_array[3,i] + op_array[15,i] * ΔT, g)
end
function calc_coord!(op_array, phytos, arch::Architecture, g::Grids, ΔT)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, g, ΔT)
    wait(device(arch), event)
    return nothing
end
