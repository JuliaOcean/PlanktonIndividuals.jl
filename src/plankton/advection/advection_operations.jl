##### deal with particles moved out of the domain
function periodic_domain!(plank, ac, g::Grids)
    plank.x .= isless.(plank.x, g.xF[g.Nx+g.Hx+1]) .* plank.x .+
               isequal.(plank.x, g.xF[g.Nx+g.Hx+1]) .* plank.x .+
               isless.(g.xF[g.Nx+g.Hx+1], plank.x) .* (plank.x .- (g.Nx*g.Δx))

    plank.y .= isless.(plank.y, g.yF[g.Ny+g.Hy+1]) .* plank.y .+
               isequal.(plank.y, g.yF[g.Ny+g.Hy+1]) .* plank.y .+
               isless.(g.yF[g.Ny+g.Hy+1], plank.y) .* (plank.y .- (g.Ny*g.Δy))

    # bounded for z direction
    plank.z .= isless.(plank.z, g.zF[g.Nz+g.Hz+1]) .* plank.z .+
               isequal.(plank.z, g.zF[g.Nz+g.Hz+1]) .* plank.z .+
               isless.(g.zF[g.Nz+g.Hz+1], plank.z) .* (plank.z .- 0.01)

    plank.x .= isless.(plank.x, g.xF[g.Hx+1]) .* (plank.x .+ (g.Nx*g.Δx)) .+
               isequal.(plank.x, g.xF[g.Hx+1]) .* plank.x .+
               isless.(g.xF[g.Hx+1], plank.x) .* plank.x

    plank.y .= isless.(plank.y, g.yF[g.Hy+1]) .* (plank.y .+ (g.Ny*g.Δy)) .+
               isequal.(plank.y, g.yF[g.Hy+1]) .* plank.y .+
               isless.(g.yF[g.Hy+1], plank.y) .* plank.y

    # bounded for z direction
    plank.z .= isless.(plank.z, g.zF[g.Hz+1]) .* plank.z .+
               isequal.(plank.z, g.zF[g.Hz+1]) .* plank.z .+
               isless.(g.zF[g.Hz+1], plank.z) .* plank.z

    plank.x .*= ac
    plank.y .*= ac
    plank.z .*= ac

    return nothing
end

##### find indices (halo points excluded)
function find_inds!(plank, ac, g::Grids)
    plank.xi .= (plank.x .- g.xF[g.Hx+1]) .÷ g.Δx .+ 1
    plank.yi .= (plank.y .- g.yF[g.Hy+1]) .÷ g.Δy .+ 1
    plank.zi .= (plank.z .- g.zF[g.Hz+1]) .÷ g.Δz .+ 1

    plank.xi .*= ac
    plank.yi .*= ac
    plank.zi .*= ac

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
