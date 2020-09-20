##### deal with particles moved out of the domain
@kernel function periodic_domain⁺_kernel!(plank, g::Grids, ind)
    i = @index(Global, Linear)
    if plank[i,1+ind] ≥ g.xF[g.Nx+g.Hx+1]
        plank[i,1+ind] = plank[i,1+ind] - g.Nx*g.Δx
        # plank[i,1+ind] = g.xF[g.Nx+g.Hx+1]
    end
    if plank[i,2+ind] ≥ g.yF[g.Ny+g.Hy+1]
        plank[i,2+ind] = plank[i,2+ind] - g.Ny*g.Δy
        # plank[i,2+ind] = g.yF[g.Ny+g.Hy+1]
    end
    if plank[i,3+ind] ≥ g.zF[g.Nz+g.Hz+1]
        plank[i,3+ind] = g.zF[g.Nz+g.Hz+1] - 1.0e-10 # to keep in the boundary
    end
end
function periodic_domain⁺!(plank, arch::Architecture, g::Grids, ind::Int64)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = periodic_domain⁺_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, g, ind)
    wait(device(arch), event)
    return nothing
end
@kernel function periodic_domain⁻_kernel!(plank, g::Grids, ind)
    i = @index(Global, Linear)
    if plank[i,1+ind] < g.xF[g.Hx+1]
        # plank[i,1+ind] = g.xF[g.Hx+1]
        plank[i,1+ind] = plank[i,1+ind] + g.Nx*g.Δx
    end
    if plank[i,2+ind] < g.yF[g.Hy+1]
        # plank[i,2+ind] = g.yF[g.Hy+1]
        plank[i,2+ind] = plank[i,2+ind] + g.Ny*g.Δy
    end
    if plank[i,3+ind] < g.zF[g.Hz+1]
        plank[i,3+ind] = g.zF[g.Hz+1]
    end
end
function periodic_domain⁻!(plank, arch::Architecture, g::Grids, ind::Int64)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = periodic_domain⁻_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, g, ind)
    wait(device(arch), event)
    return nothing
end
function in_domain!(plank, arch::Architecture, g::Grids, ind)
    periodic_domain⁺!(plank, arch::Architecture, g::Grids, ind)
    periodic_domain⁻!(plank, arch::Architecture, g::Grids, ind)
end

in_domain!(plank, arch::Architecture, g::Grids) = in_domain!(plank, arch::Architecture, g::Grids, 0)


##### find indices (halo points excluded)
@inline find_xF_ind(x, g::Grids) = x ≥ g.xF[g.Hx+1] ? x ÷ g.Δx + 1 : 0.0

@inline find_yF_ind(y, g::Grids) = y ≥ g.yF[g.Hy+1] ? y ÷ g.Δy + 1 : 0.0

@inline find_zF_ind(z, g::Grids) = z ≥ g.zF[g.Hz+1] ? (g.Nz * g.Δz + z) ÷ g.Δz + 1 : 0.0

##### find cell and face indices for each individual
@kernel function find_inds_kernel!(plank, g::Grids, indₜ::Int64, ind₀::Int64)
    i = @index(Global, Linear)
    @inbounds plank[i,1+indₜ] = find_xF_ind(plank[i,1+ind₀], g)  # xF index
    @inbounds plank[i,2+indₜ] = find_yF_ind(plank[i,2+ind₀], g)  # yF index
    @inbounds plank[i,3+indₜ] = find_zF_ind(plank[i,3+ind₀], g)  # zF index
end
function find_inds!(plank, arch::Architecture, g::Grids, indₜ::Int64, ind₀::Int64)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = find_inds_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, g, indₜ, ind₀)
    wait(device(arch), event)
    return nothing
end

##### find velocities around the individual to interpolate
@kernel function find_vel_kernel!(plank, ind_array::AbstractArray{Int64,2}, g::Grids, u, v, w)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds plank[i,37] = u[x₀,   y₀,   z₀  ]
    @inbounds plank[i,38] = u[x₀+1, y₀,   z₀  ]
    @inbounds plank[i,39] = v[x₀,   y₀,   z₀  ]
    @inbounds plank[i,40] = v[x₀,   y₀+1, z₀  ]
    @inbounds plank[i,41] = w[x₀,   y₀,   z₀  ]
    @inbounds plank[i,42] = w[x₀,   y₀,   z₀+1]
end
function find_vel!(plank, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids, u, v, w)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = find_vel_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, ind_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

@kernel function find_xᵈ_kernel!(plank, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    x₀ = ind_array[i,1] + g.Hx
    y₀ = ind_array[i,2] + g.Hy
    z₀ = ind_array[i,3] + g.Hz
    @inbounds plank[i,43] = (plank[i,1] - g.xF[x₀]) / g.Δx
    @inbounds plank[i,44] = (plank[i,2] - g.yF[y₀]) / g.Δy
    @inbounds plank[i,45] = (plank[i,3] - g.zF[z₀]) / g.Δz
end
function find_xᵈ!(plank, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = find_xᵈ_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, ind_array, g)
    wait(device(arch), event)
    return nothing
end

##### velocity interpolation for each individual
@inline linear_itpl(u0, u1, xd) = u0 * (1.0 - xd) + u1 * xd

@kernel function vel_interpolation_kernel!(plank, g::Grids, ind::Int64)
    i = @index(Global, Linear)
    @inbounds plank[i, 46+ind*3] = linear_itpl(plank[i,37], plank[i,38], plank[i,43])
    @inbounds plank[i, 47+ind*3] = linear_itpl(plank[i,39], plank[i,40], plank[i,44])
    @inbounds plank[i, 48+ind*3] = linear_itpl(plank[i,41], plank[i,42], plank[i,45])
end
function vel_interpolation!(plank, arch::Architecture, g::Grids, ind::Int64)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = vel_interpolation_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, g, ind)
    wait(device(arch), event)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_1st_intermediate_coord_kernel!(plank, ΔT)
    i = @index(Global, Linear)
    @inbounds plank[i,34] = plank[i,1] + 0.5 * plank[i,46] * ΔT
    @inbounds plank[i,35] = plank[i,2] + 0.5 * plank[i,47] * ΔT
    @inbounds plank[i,36] = plank[i,3] + 0.5 * plank[i,48] * ΔT
end
@kernel function calc_2nd_intermediate_coord_kernel!(plank, ΔT)
    i = @index(Global, Linear)
    @inbounds plank[i,34] = plank[i,34] + 0.5 * plank[i,49] * ΔT
    @inbounds plank[i,35] = plank[i,35] + 0.5 * plank[i,50] * ΔT
    @inbounds plank[i,36] = plank[i,36] + 0.5 * plank[i,51] * ΔT
end
@kernel function calc_3rd_intermediate_coord_kernel!(plank, ΔT)
    i = @index(Global, Linear)
    @inbounds plank[i,34] = plank[i,34] + 1.0 * plank[i,52] * ΔT
    @inbounds plank[i,35] = plank[i,35] + 1.0 * plank[i,53] * ΔT
    @inbounds plank[i,36] = plank[i,36] + 1.0 * plank[i,54] * ΔT
end
function calc_1st_intermediate_coord!(plank, arch::Architecture, ΔT)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = calc_1st_intermediate_coord_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, ΔT)
    wait(device(arch), event)
    return nothing
end
function calc_2nd_intermediate_coord!(plank, arch::Architecture, ΔT)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = calc_2nd_intermediate_coord_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, ΔT)
    wait(device(arch), event)
    return nothing
end
function calc_3rd_intermediate_coord!(plank, arch::Architecture, ΔT)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = calc_3rd_intermediate_coord_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate final velocities by RK4
@kernel function calc_vel_rk4_kernel!(plank)
    i = @index(Global, Linear)
    @inbounds plank[i,46] = (plank[i,46] + 2*plank[i,49] + 2*plank[i,52] + plank[i,55]) /6
    @inbounds plank[i,47] = (plank[i,47] + 2*plank[i,50] + 2*plank[i,53] + plank[i,56]) /6
    @inbounds plank[i,48] = (plank[i,48] + 2*plank[i,51] + 2*plank[i,54] + plank[i,57]) /6
end
function calc_vel_rk4!(plank, arch::Architecture)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = calc_vel_rk4_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank)
    wait(device(arch), event)
    return nothing
end

##### calculate coordinates of each individual
@kernel function calc_coord_kernel!(plank, ΔT)
    i = @index(Global, Linear)
    @inbounds plank[i,1] = plank[i,1] + plank[i,46] * ΔT
    @inbounds plank[i,2] = plank[i,2] + plank[i,47] * ΔT
    @inbounds plank[i,3] = plank[i,3] + plank[i,48] * ΔT
end
function calc_coord!(plank, arch::Architecture, ΔT)
    plank_num = floor(Int64, sum(plank[:,61]))
    kernel! = calc_coord_kernel!(device(arch), 256, (plank_num,))
    event = kernel!(plank, ΔT)
    wait(device(arch), event)
    return nothing
end
