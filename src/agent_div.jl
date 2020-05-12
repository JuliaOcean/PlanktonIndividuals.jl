#########################################################################
# advection  of agents (use Interpolations.jl to interpolate velocity fields)#
#########################################################################
"""
    generate_vel_itp(grid, vel)
Use Interpolations.jl to generate interpolation objects
"""
function generate_vel_itp(grid, vel)
    u_itp = interpolate((grid.xF,grid.yC,grid.zC),vel.u,Gridded(Linear()))
    v_itp = interpolate((grid.xC,grid.yF,grid.zC),vel.v,Gridded(Linear()))
    w_itp = interpolate((grid.xC,grid.yC,grid.zF),vel.w,Gridded(Linear()))
    return (u_itp, v_itp, w_itp)
end

"""
    get_vels(x, y, z, vel_itp)
Read and interpolate velocity at '(x,y,z)'
'(x,y,z)' is actual location of an individual
'vel_itp' is interpolation objects
"""
function get_vels(x, y, z, vel_itp)
    u_itp = vel_itp[1]; v_itp = vel_itp[2]; w_itp = vel_itp[3];
    u = u_itp(x, y, z); v = v_itp(x, y, z); w = w_itp(x, y, z);
    return u, v, w
end

"""
    periodic_domain(xF, x)
'xF' is the array containing the faces of grids
'x' is the coordinate of an individual
"""
function periodic_domain(xF, x)
    xF₀ = xF[2]; xFₜ= xF[end]
    if xF₀ < x < xFₜ
        return x
    elseif x ≥ xFₜ
        return x - xFₜ
    elseif x ≤ xF₀
        return x + xFₜ
    end
end

"""
    agent_advection(phyts_a,vel_itp,g,ΔT)
Update grid indices of all the individuals according to velocity fields of each time step
Periodic domain is used
'phyts_a' is a dataframe contains all the individuals of current time step
'vel_itp' is the tuple contains interpolations of u, v, w velocites of current time step
'g' is the grid information and 'ΔT' is time step
"""
function agent_advection(p_xyz,vel_itp,g,ΔT::Int64)
    uvel, vvel, wvel = get_vels(p_xyz[1], p_xyz[2], p_xyz[3], vel_itp)
    p_xyz[1] = p_xyz[1] + uvel*ΔT
    p_xyz[2] = p_xyz[2] + vvel*ΔT
    p_xyz[3] = max(g.zF[2],min(g.zF[end-1],p_xyz[3] + wvel*ΔT))
    # periodic domain
    p_xyz[1] = periodic_domain(g.xF, p_xyz[1])
    p_xyz[2] = periodic_domain(g.yF, p_xyz[2])
    return p_xyz
end
"""
    agent_advectionRK4(p_xyzs_a, vel_itps, g, ΔT::Int64)
'vel_itps' is an array of tuples containing interpolations of u, v, w velocites of current time step
"""
function agent_advectionRK4(p_xyz, vel_itps, g, ΔT::Int64)
    u1,v1,w1 = get_vels(p_xyz[1], p_xyz[2], p_xyz[3], vel_itps[1]) # velocites at t
    gx1 = periodic_domain(g.xF, p_xyz[1] + u1*0.5*ΔT)
    gy1 = periodic_domain(g.yF, p_xyz[2] + v1*0.5*ΔT)
    gz1 = p_xyz[3] + w1*0.5*ΔT
    gz1 = max(g.zF[2],min(g.zF[end-1],gz1))
    u2,v2,w2 = get_vels(gx1, gy1, gz1, vel_itps[2]) # velocites at t+0.5ΔT
    gx2 = periodic_domain(g.xF, p_xyz[1] + u2*0.5*ΔT)
    gy2 = periodic_domain(g.yF, p_xyz[2] + v2*0.5*ΔT)
    gz2 = p_xyz[3] + w2*0.5*ΔT
    gz2 = max(g.zF[2],min(g.zF[end-1],gz2))
    u3,v3,w3 = get_vels(gx2, gy2, gz2, vel_itps[2]) # velocites at t+0.5ΔT
    gx3 = periodic_domain(g.xF, p_xyz[1] + u3*0.5*ΔT)
    gy3 = periodic_domain(g.yF, p_xyz[2] + v3*0.5*ΔT)
    gz3 = p_xyz[3] + w3*0.5*ΔT
    gz3 = max(g.zF[2],min(g.zF[end-1],gz3))
    u4,v4,w4 = get_vels(gx3, gy3, gz3, vel_itps[3]) # velocites at t+ΔT
    dx = (u1 + 2*u2 + 2*u3 + u4) / 6 * ΔT
    dy = (v1 + 2*v2 + 2*v3 + v4) / 6 * ΔT
    dz = (w1 + 2*w2 + 2*w3 + w4) / 6 * ΔT
    p_xyz[1] = periodic_domain(g.xF, p_xyz[1] + dx)
    p_xyz[2] = periodic_domain(g.yF, p_xyz[2] + dy)
    p_xyz[3] = p_xyz[3] + dz
    p_xyz[3] = max(g.zF[2], min(g.zF[end-1], p_xyz[3]))
    return p_xyz
end

"""
    agent_diffusionX(phyt,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionX(p_x, g, κh, ΔT::Int64)
    p_x += rand(Uniform(-1.0,1.0)) * κh * ΔT
    p_x = periodic_domain(g.xF, p_x)
    return p_x
end

"""
    agent_diffusionY(phyt,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionY(p_y, g, κh, ΔT::Int64)
    p_y += rand(Uniform(-1.0,1.0)) * κh * ΔT
    p_y = periodic_domain(g.yF, p_y)
    return p_y
end
"""
    agent_diffusionZ(phyt,g,κv)
Using a random walk algorithm for vertical diffusion
"""
function agent_diffusionZ(p_z, g, κv, ΔT::Int64)
    p_z += rand(Uniform(-1.0,1.0)) * κv * ΔT
    p_z = max(g.zF[2], min(g.zF[end-1], p_z))
    return p_z
end
