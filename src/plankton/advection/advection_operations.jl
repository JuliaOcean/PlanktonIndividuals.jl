##### deal with particles moved out of the domain
function periodic_domain!(plank, ac, g::Grids)
    plank.x .= isless.(plank.x, g.xF[g.x⁺]) .* plank.x .+
               isequal.(plank.x, g.xF[g.x⁺]) .* plank.x .+
               isless.(g.xF[g.x⁺], plank.x) .* (plank.x .- (g.Nx*g.Δx))

    plank.y .= isless.(plank.y, g.yF[g.y⁺]) .* plank.y .+
               isequal.(plank.y, g.yF[g.y⁺]) .* plank.y .+
               isless.(g.yF[g.y⁺], plank.y) .* (plank.y .- (g.Ny*g.Δy))

    plank.z .= isless.(plank.z, g.zF[g.z⁺]) .* plank.z .+
               isequal.(plank.z, g.zF[g.z⁺]) .* plank.z .+
               isless.(g.zF[g.z⁺], plank.z) .* (plank.z .- (g.Nz*g.Δz))

    plank.x .= isless.(plank.x, g.xF[g.x⁻]) .* (plank.x .+ (g.Nx*g.Δx)) .+
               isequal.(plank.x, g.xF[g.x⁻]) .* plank.x .+
               isless.(g.xF[g.x⁻], plank.x) .* plank.x

    plank.y .= isless.(plank.y, g.yF[g.y⁻]) .* (plank.y .+ (g.Ny*g.Δy)) .+
               isequal.(plank.y, g.yF[g.y⁻]) .* plank.y .+
               isless.(g.yF[g.y⁻], plank.y) .* plank.y

    plank.z .= isless.(plank.z, g.zF[g.z⁻]) .* (plank.z .+ (g.Nz*g.Δz)) .+
               isequal.(plank.z, g.zF[g.z⁻]) .* plank.z .+
               isless.(g.zF[g.z⁻], plank.z) .* plank.z

    plank.x .*= ac
    plank.y .*= ac
    plank.z .*= ac

    return nothing
end

##### find indices (halo points excluded)
# indices of x and y must start from 0
# indices of z must end at 0
function find_inds!(plank, coord, ac, g::Grids)
    coord.x .= plank.x .÷ g.Δx .+ 1
    coord.y .= plank.y .÷ g.Δy .+ 1
    coord.z .= (plank.z .+ (g.Nz*g.Δz)) .÷ g.Δz .+ 1

    coord.x .*= ac
    coord.y .*= ac
    coord.z .*= ac

    return nothing
end

##### find velocities around the individual to interpolate
function find_vel!(velos, x, y, z, ac, g::Grids, u, v, w)
    @inbounds velos.u⁻ .= u[CartesianIndex.(x .+ g.Hx,      y .+ g.Hy,      z .+ g.Hz     )] .* ac
    @inbounds velos.u⁺ .= u[CartesianIndex.(x .+ g.Hx .+ 1, y .+ g.Hy,      z .+ g.Hz     )] .* ac
    @inbounds velos.v⁻ .= v[CartesianIndex.(x .+ g.Hx,      y .+ g.Hy,      z .+ g.Hz     )] .* ac
    @inbounds velos.v⁺ .= v[CartesianIndex.(x .+ g.Hx,      y .+ g.Hy .+ 1, z .+ g.Hz     )] .* ac
    @inbounds velos.w⁻ .= w[CartesianIndex.(x .+ g.Hx,      y .+ g.Hy,      z .+ g.Hz     )] .* ac
    @inbounds velos.w⁺ .= w[CartesianIndex.(x .+ g.Hx,      y .+ g.Hy,      z .+ g.Hz .+ 1)] .* ac
    return nothing
end

function find_xᵈ!(velos, plank, g::Grids)
    velos.xd .= plank.x .% g.Δx ./ g.Δx
    velos.yd .= plank.y .% g.Δy ./ g.Δy
    velos.zd .= plank.z .% g.Δz ./ g.Δz
    return nothing
end

##### velocity interpolation for each individual
@inline linear_itpl(u0, u1, xd) = u0 * (1.0 - xd) + u1 * xd

function vel1_interpolation!(velos)
    velos.u1 .= linear_itpl.(velos.u⁻, velos.u⁺, velos.xd)
    velos.v1 .= linear_itpl.(velos.v⁻, velos.v⁺, velos.yd)
    velos.w1 .= linear_itpl.(velos.w⁻, velos.w⁺, velos.zd)
    return nothing
end
function vel2_interpolation!(velos)
    velos.u2 .= linear_itpl.(velos.u⁻, velos.u⁺, velos.xd)
    velos.v2 .= linear_itpl.(velos.v⁻, velos.v⁺, velos.yd)
    velos.w2 .= linear_itpl.(velos.w⁻, velos.w⁺, velos.zd)
    return nothing
end
function vel3_interpolation!(velos)
    velos.u3 .= linear_itpl.(velos.u⁻, velos.u⁺, velos.xd)
    velos.v3 .= linear_itpl.(velos.v⁻, velos.v⁺, velos.yd)
    velos.w3 .= linear_itpl.(velos.w⁻, velos.w⁺, velos.zd)
    return nothing
end
function vel4_interpolation!(velos)
    velos.u4 .= linear_itpl.(velos.u⁻, velos.u⁺, velos.xd)
    velos.v4 .= linear_itpl.(velos.v⁻, velos.v⁺, velos.yd)
    velos.w4 .= linear_itpl.(velos.w⁻, velos.w⁺, velos.zd)
    return nothing
end

##### calculate intermediate coordinates
function calc_coord_1!(plank, velos, ΔT)
    velos.x .= plank.x .+ 0.5 .* velos.u1 .* ΔT .* plank.ac
    velos.y .= plank.y .+ 0.5 .* velos.v1 .* ΔT .* plank.ac
    velos.z .= plank.z .+ 0.5 .* velos.w1 .* ΔT .* plank.ac
    return nothing
end
function calc_coord_2!(plank, velos, ΔT)
    velos.x .= plank.x .+ 0.5 .* velos.u2 .* ΔT .* plank.ac
    velos.y .= plank.y .+ 0.5 .* velos.v2 .* ΔT .* plank.ac
    velos.z .= plank.z .+ 0.5 .* velos.w2 .* ΔT .* plank.ac
    return nothing
end
function calc_coord_3!(plank, velos, ΔT)
    velos.x .= plank.x .+ 1.0 .* velos.u3 .* ΔT .* plank.ac
    velos.y .= plank.y .+ 1.0 .* velos.v3 .* ΔT .* plank.ac
    velos.z .= plank.z .+ 1.0 .* velos.w3 .* ΔT .* plank.ac
    return nothing
end

##### calculate final velocities by RK4
function calc_vel_rk4!(velos)
    velos.u1 .= (velos.u1 .+ 2 .* velos.u2 .+ 2 .* velos.u3 .+ velos.u4) ./ 6
    velos.v1 .= (velos.v1 .+ 2 .* velos.v2 .+ 2 .* velos.v3 .+ velos.v4) ./ 6
    velos.w1 .= (velos.w1 .+ 2 .* velos.w2 .+ 2 .* velos.w3 .+ velos.w4) ./ 6
    return nothing
end

##### calculate coordinates of each individual
function calc_coord!(plank, velos, ΔT)
    plank.x .= plank.x .+ velos.u1 .* ΔT .* plank.ac
    plank.y .= plank.y .+ velos.v1 .* ΔT .* plank.ac
    plank.z .= plank.z .+ velos.w1 .* ΔT .* plank.ac
    return nothing
end

# @kernel function vel_interpolation_kernel!(plank, ind::Int64)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         @inbounds plank[i, 46+ind*3] = linear_itpl(plank[i,37], plank[i,38], plank[i,43])
#         @inbounds plank[i, 47+ind*3] = linear_itpl(plank[i,39], plank[i,40], plank[i,44])
#         @inbounds plank[i, 48+ind*3] = linear_itpl(plank[i,41], plank[i,42], plank[i,45])
#     end
# end
# function vel_interpolation!(plank, arch::Architecture, ind::Int64)
#     kernel! = vel_interpolation_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, ind)
#     wait(device(arch), event)
#     return nothing
# end

# ##### calculate intermediate coordinates
# @kernel function calc_intermediate_coord_kernel!(plank, ΔT, weight::Float64, ind::Int64)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         @inbounds plank[i,34] = plank[i,1] + weight * plank[i,46+ind*3] * ΔT
#         @inbounds plank[i,35] = plank[i,2] + weight * plank[i,47+ind*3] * ΔT
#         @inbounds plank[i,36] = plank[i,3] + weight * plank[i,48+ind*3] * ΔT
#     end
# end
# function calc_intermediate_coord!(plank, arch::Architecture, ΔT, weight::Float64, ind::Int64)
#     kernel! = calc_intermediate_coord_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, ΔT, weight, ind)
#     wait(device(arch), event)
#     return nothing
# end

# ##### calculate final velocities by RK4
# @kernel function calc_vel_rk4_kernel!(plank)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         @inbounds plank[i,46] = (plank[i,46] + 2*plank[i,49] + 2*plank[i,52] + plank[i,55]) /6
#         @inbounds plank[i,47] = (plank[i,47] + 2*plank[i,50] + 2*plank[i,53] + plank[i,56]) /6
#         @inbounds plank[i,48] = (plank[i,48] + 2*plank[i,51] + 2*plank[i,54] + plank[i,57]) /6
#     end
# end
# function calc_vel_rk4!(plank, arch::Architecture)
#     kernel! = calc_vel_rk4_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank)
#     wait(device(arch), event)
#     return nothing
# end

# ##### calculate coordinates of each individual
# @kernel function calc_coord_kernel!(plank, ΔT)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         @inbounds plank[i,1] = plank[i,1] + plank[i,46] * ΔT
#         @inbounds plank[i,2] = plank[i,2] + plank[i,47] * ΔT
#         @inbounds plank[i,3] = plank[i,3] + plank[i,48] * ΔT
#     end
# end
# function calc_coord!(plank, arch::Architecture, ΔT)
#     kernel! = calc_coord_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, ΔT)
#     wait(device(arch), event)
#     return nothing
# end

# @kernel function periodic_domain⁺_kernel!(plank, g::Grids, ind)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         if plank[i,1+ind] > g.xF[g.Nx+g.Hx+1]
#             @inbounds plank[i,1+ind] = plank[i,1+ind] - g.Nx*g.Δx
#         end
#         if plank[i,2+ind] > g.yF[g.Ny+g.Hy+1]
#             @inbounds plank[i,2+ind] = plank[i,2+ind] - g.Ny*g.Δy
#         end
#         if plank[i,3+ind] ≥ g.zF[g.Nz+g.Hz+1]
#             @inbounds plank[i,3+ind] = g.zF[g.Nz+g.Hz+1] - 1.0e-2 # to keep in the boundary
#         end
#     end
# end
# function periodic_domain⁺!(plank, arch::Architecture, g::Grids, ind::Int64)
#     kernel! = periodic_domain⁺_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, g, ind)
#     wait(device(arch), event)
#     return nothing
# end
# @kernel function periodic_domain⁻_kernel!(plank, g::Grids, ind)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         if plank[i,1+ind] < g.xF[g.Hx+1]
#             @inbounds plank[i,1+ind] = plank[i,1+ind] + g.Nx*g.Δx
#         end
#         if plank[i,2+ind] < g.yF[g.Hy+1]
#             @inbounds plank[i,2+ind] = plank[i,2+ind] + g.Ny*g.Δy
#         end
#         if plank[i,3+ind] < g.zF[g.Hz+1]
#             @inbounds plank[i,3+ind] = g.zF[g.Hz+1]
#         end
#     end
# end
# function periodic_domain⁻!(plank, arch::Architecture, g::Grids, ind::Int64)
#     kernel! = periodic_domain⁻_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, g, ind)
#     wait(device(arch), event)
#     return nothing
# end
# function in_domain!(plank, arch::Architecture, g::Grids, ind)
#     periodic_domain⁺!(plank, arch::Architecture, g::Grids, ind)
#     periodic_domain⁻!(plank, arch::Architecture, g::Grids, ind)
# end

# in_domain!(plank, arch::Architecture, g::Grids) = in_domain!(plank, arch::Architecture, g::Grids, 0)

# @inline find_xF_ind(x, g::Grids) = x ≥ g.xF[g.Hx+1] ? x ÷ g.Δx + 1 : 0.0

# @inline find_yF_ind(y, g::Grids) = y ≥ g.yF[g.Hy+1] ? y ÷ g.Δy + 1 : 0.0

# @inline find_zF_ind(z, g::Grids) = z ≥ g.zF[g.Hz+1] ? (g.Nz * g.Δz + z) ÷ g.Δz + 1 : 0.0

# ##### find cell and face indices for each individual
# @kernel function find_inds_kernel!(plank, g::Grids, ind::Int64)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         @inbounds plank[i,13] = find_xF_ind(plank[i,1+ind], g)  # xF index
#         @inbounds plank[i,14] = find_yF_ind(plank[i,2+ind], g)  # yF index
#         @inbounds plank[i,15] = find_zF_ind(plank[i,3+ind], g)  # zF index
#     end
# end
# function find_inds!(plank, arch::Architecture, g::Grids, ind::Int64)
#     kernel! = find_inds_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, g, ind)
#     wait(device(arch), event)
#     return nothing
# end

# @kernel function find_xᵈ_kernel!(velos, x, y, z, ac, plank, g::Grids)
#     i = @index(Global, Linear)
#     if ac[i] == 1.0
#         @inbounds x₀ = x[i] + g.Hx
#         @inbounds y₀ = y[i] + g.Hy
#         @inbounds z₀ = z[i] + g.Hz
#         @inbounds velos.xd[i] = (plank.x[i] - g.xF[x₀]) / g.Δx
#         @inbounds velos.yd[i] = (plank.y[i] - g.yF[y₀]) / g.Δy
#         @inbounds velos.zd[i] = (plank.z[i] - g.zF[z₀]) / g.Δz
#     end
# end
# function find_xᵈ!(velos, x, y, z, ac, plank, arch::Architecture, g::Grids)
#     kernel! = find_xᵈ_kernel!(device(arch), 256, (size(ac,1),))
#     event = kernel!(velos, x, y, z, ac, plank, g)
#     wait(device(arch), event)
#     return nothing
# end

# @kernel function find_vel_kernel!(velos, coord, x, y, z, ac, g::Grids, u, v, w)
#     i = @index(Global, Linear)
#     if ac[i] == 1.0
#         @inbounds x₀ = x[i] + g.Hx
#         @inbounds y₀ = y[i] + g.Hy
#         @inbounds z₀ = z[i] + g.Hz
#         @inbounds velos.u⁻[i] = u[x₀,   y₀,   z₀  ]
#         @inbounds velos.u⁺[i] = u[x₀+1, y₀,   z₀  ]
#         @inbounds velos.v⁻[i] = v[x₀,   y₀,   z₀  ]
#         @inbounds velos.v⁺[i] = v[x₀,   y₀+1, z₀  ]
#         @inbounds velos.w⁻[i] = w[x₀,   y₀,   z₀  ]
#         @inbounds velos.w⁺[i] = w[x₀,   y₀,   z₀+1]
#         @inbounds coord.xF[i] = g.xF[x₀]
#         @inbounds coord.yF[i] = g.yF[y₀]
#         @inbounds coord.zF[i] = g.zF[z₀]
#     end
# end
# function find_vel!(velos, coord, x, y, z, ac, arch::Architecture, g::Grids, u, v, w)
#     kernel! = find_vel_kernel!(device(arch), 256, (size(ac,1),))
#     event = kernel!(velos, coord, x, y, z, ac, g, u, v, w)
#     wait(device(arch), event)
#     return nothing
# end
