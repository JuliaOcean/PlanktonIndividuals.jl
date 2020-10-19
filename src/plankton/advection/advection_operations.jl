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

