##### deal with particles moved out of the domain
@kernel function periodic_domain⁺_kernel!(plank, g::Grids, ind)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,1+ind] > g.xF[g.Nx+g.Hx+1]
            plank[i,1+ind] = plank[i,1+ind] - g.Nx*g.Δx
        end
        if plank[i,2+ind] > g.yF[g.Ny+g.Hy+1]
            plank[i,2+ind] = plank[i,2+ind] - g.Ny*g.Δy
        end
        if plank[i,3+ind] ≥ g.zF[g.Nz+g.Hz+1]
            plank[i,3+ind] = g.zF[g.Nz+g.Hz+1] - 1.0e-2 # to keep in the boundary
        end
    end
end
function periodic_domain⁺!(plank, arch::Architecture, g::Grids, ind::Int64)
    kernel! = periodic_domain⁺_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, g, ind)
    wait(device(arch), event)
    return nothing
end
@kernel function periodic_domain⁻_kernel!(plank, g::Grids, ind)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,1+ind] < g.xF[g.Hx+1]
            plank[i,1+ind] = plank[i,1+ind] + g.Nx*g.Δx
        end
        if plank[i,2+ind] < g.yF[g.Hy+1]
            plank[i,2+ind] = plank[i,2+ind] + g.Ny*g.Δy
        end
        if plank[i,3+ind] < g.zF[g.Hz+1]
            plank[i,3+ind] = g.zF[g.Hz+1]
        end
    end
end
function periodic_domain⁻!(plank, arch::Architecture, g::Grids, ind::Int64)
    kernel! = periodic_domain⁻_kernel!(device(arch), 256, (size(plank,1),))
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
    if plank[i,58] == 1.0
        @inbounds plank[i,1+indₜ] = find_xF_ind(plank[i,1+ind₀], g)  # xF index
        @inbounds plank[i,2+indₜ] = find_yF_ind(plank[i,2+ind₀], g)  # yF index
        @inbounds plank[i,3+indₜ] = find_zF_ind(plank[i,3+ind₀], g)  # zF index
    end
end
function find_inds!(plank, arch::Architecture, g::Grids, indₜ::Int64, ind₀::Int64)
    kernel! = find_inds_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, g, indₜ, ind₀)
    wait(device(arch), event)
    return nothing
end

##### find velocities around the individual to interpolate
@kernel function find_vel_kernel!(plank, ind_array::AbstractArray{Int64,2}, g::Grids, u, v, w)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds x₀ = ind_array[i,1] + g.Hx
        @inbounds y₀ = ind_array[i,2] + g.Hy
        @inbounds z₀ = ind_array[i,3] + g.Hz
        @inbounds plank[i,37] = u[x₀,   y₀,   z₀  ]
        @inbounds plank[i,38] = u[x₀+1, y₀,   z₀  ]
        @inbounds plank[i,39] = v[x₀,   y₀,   z₀  ]
        @inbounds plank[i,40] = v[x₀,   y₀+1, z₀  ]
        @inbounds plank[i,41] = w[x₀,   y₀,   z₀  ]
        @inbounds plank[i,42] = w[x₀,   y₀,   z₀+1]
    end
end
function find_vel!(plank, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids, u, v, w)
    kernel! = find_vel_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, ind_array, g, u, v, w)
    wait(device(arch), event)
    return nothing
end

@kernel function find_xᵈ_kernel!(plank, ind_array::AbstractArray{Int64,2}, g::Grids)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds x₀ = ind_array[i,1] + g.Hx
        @inbounds y₀ = ind_array[i,2] + g.Hy
        @inbounds z₀ = ind_array[i,3] + g.Hz
        @inbounds plank[i,43] = (plank[i,1] - g.xF[x₀]) / g.Δx
        @inbounds plank[i,44] = (plank[i,2] - g.yF[y₀]) / g.Δy
        @inbounds plank[i,45] = (plank[i,3] - g.zF[z₀]) / g.Δz
    end
end
function find_xᵈ!(plank, ind_array::AbstractArray{Int64,2}, arch::Architecture, g::Grids)
    kernel! = find_xᵈ_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, ind_array, g)
    wait(device(arch), event)
    return nothing
end

##### velocity interpolation for each individual
@inline linear_itpl(u0, u1, xd) = u0 * (1.0 - xd) + u1 * xd

function vel_interpolation!(plank, ind::Int64, active_num::Int64)
    plank[1:active_num, 46+ind*3] .= linear_itpl.(plank[1:active_num,37], plank[1:active_num,38], plank[1:active_num,43])
    plank[1:active_num, 47+ind*3] .= linear_itpl.(plank[1:active_num,39], plank[1:active_num,40], plank[1:active_num,44])
    plank[1:active_num, 48+ind*3] .= linear_itpl.(plank[1:active_num,41], plank[1:active_num,42], plank[1:active_num,45])
end


##### calculate intermediate coordinates
function calc_intermediate_coord!(plank, ΔT, weight::Float64, active_num::Int64)
    plank[1:active_num,34] .= plank[1:active_num,1] .+ weight .* plank[1:active_num,46] .* ΔT
    plank[1:active_num,35] .= plank[1:active_num,2] .+ weight .* plank[1:active_num,47] .* ΔT
    plank[1:active_num,36] .= plank[1:active_num,3] .+ weight .* plank[1:active_num,48] .* ΔT
end

##### calculate final velocities by RK4
function calc_vel_rk4!(plank, active_num::Int64)
    plank[1:active_num,46] .= (plank[1:active_num,46] .+ 2 .* plank[1:active_num,49] .+
                               2 .* plank[1:active_num,52] .+ plank[1:active_num,55]) ./ 6
    plank[1:active_num,47] .= (plank[1:active_num,47] .+ 2 .* plank[1:active_num,50] .+
                               2 .* plank[1:active_num,53] .+ plank[1:active_num,56]) ./ 6
    plank[1:active_num,48] .= (plank[1:active_num,48] .+ 2 .* plank[1:active_num,51] .+
                               2 .* plank[1:active_num,54] .+ plank[1:active_num,57]) ./ 6
end

##### calculate coordinates of each individual
function calc_coord!(plank, ΔT, active_num::Int64)
    plank[1:active_num,1] .= plank[1:active_num,1] .+ plank[1:active_num,46] .* ΔT
    plank[1:active_num,2] .= plank[1:active_num,2] .+ plank[1:active_num,47] .* ΔT
    plank[1:active_num,3] .= plank[1:active_num,3] .+ plank[1:active_num,48] .* ΔT
end
