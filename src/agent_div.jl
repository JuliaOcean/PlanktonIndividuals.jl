#########################################################################
# advection  of agents (use Interpolations.jl to interpolate velocity fields)#
#########################################################################
"""
    generate_vel_itp(grid, vel)
Use Interpolations.jl to generate interpolation objects
"""
function generate_vel_itp(grid, vel)
    # deal with halo points
    xF = copy(grid.xF[:,1])
    pushfirst!(xF,xF[1]-(xF[2]-xF[1]))
    yF = copy(grid.yF[1,:])
    pushfirst!(yF,yF[1]-(yF[2]-yF[1]))
    zF = copy(grid.zF)
    pushfirst!(zF,zF[1]-(zF[2]-zF[1]))
    pushfirst!(zF,zF[1]-(zF[2]-zF[1]))

    xC = copy(grid.xC[:,1])
    pushfirst!(xC,xF[2]-(xC[1]-xF[2]))
    push!(xC,xF[end]+(xF[end]-xC[end]))
    yC = copy(grid.yC[1,:])
    pushfirst!(yC,yF[2]-(yC[1]-yF[2]))
    push!(yC,yF[end]+(yF[end]-yC[end]))
    zC = copy(grid.zC)
    pushfirst!(zC,zF[2]-(zC[1]-zF[2]))
    push!(zC,zF[end]+(zF[end]-zC[end]))

    u_itp = interpolate((xF,yC,zC),vel.u,Gridded(Linear()))
    v_itp = interpolate((xC,yF,zC),vel.v,Gridded(Linear()))
    w_itp = interpolate((xC,yC,zF),vel.w,Gridded(Linear()))
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
    xF₀ = xF[1]; xFₜ= xF[end]
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
function agent_advection(phyt,vel_itp,g,ΔT::Int64)
    uvel, vvel, wvel = get_vels(phyt[1], phyt[2], phyt[3], vel_itp)
    phyt[1] = phyt[1] + uvel*ΔT
    phyt[2] = phyt[2] + vvel*ΔT
    phyt[3] = max(g.zF[1],min(g.zF[end],phyt[3] + wvel*ΔT))
    # periodic domain
    phyt[1] = periodic_domain(g.xF, phyt[1])
    phyt[2] = periodic_domain(g.yF, phyt[2])
end
"""
    agent_advectionRK4(phyts_a, vel_itps, g, ΔT::Int64)
'vel_itps' is an array of tuples containing interpolations of u, v, w velocites of current time step
"""
function agent_advectionRK4(phyt, vel_itps, g, ΔT::Int64)
    u1,v1,w1 = get_vels(phyt[1], phyt[2], phyt[3], vel_itps[1]) # velocites at t
    gx1 = periodic_domain(g.xF, phyt[1] + u1*0.5*ΔT)
    gy1 = periodic_domain(g.yF, phyt[2] + v1*0.5*ΔT)
    gz1 = phyt[3] + w1*0.5*ΔT
    gz1 = max(g.zF[1],min(g.zF[end],gz1))
    u2,v2,w2 = get_vels(gx1, gy1, gz1, vel_itps[2]) # velocites at t+0.5ΔT
    gx2 = periodic_domain(g.xF, phyt[1] + u2*0.5*ΔT)
    gy2 = periodic_domain(g.yF, phyt[2] + v2*0.5*ΔT)
    gz2 = phyt[3] + w2*0.5*ΔT
    gz2 = max(g.zF[1],min(g.zF[end],gz2))
    u3,v3,w3 = get_vels(gx2, gy2, gz2, vel_itps[2]) # velocites at t+0.5ΔT
    gx3 = periodic_domain(g.xF, phyt[1] + u3*0.5*ΔT)
    gy3 = periodic_domain(g.yF, phyt[2] + v3*0.5*ΔT)
    gz3 = phyt[3] + w3*0.5*ΔT
    gz3 = max(g.zF[1],min(g.zF[end],gz3))
    u4,v4,w4 = get_vels(gx3, gy3, gz3, vel_itps[3]) # velocites at t+ΔT
    dx = (u1 + 2*u2 + 2*u3 + u4) / 6 * ΔT
    dy = (v1 + 2*v2 + 2*v3 + v4) / 6 * ΔT
    dz = (w1 + 2*w2 + 2*w3 + w4) / 6 * ΔT
    phyt[1] = periodic_domain(g.xF, phyt[1] + dx)
    phyt[2] = periodic_domain(g.yF, phyt[2] + dy)
    phyt[3] = phyt[3] + dz
    phyt[3] = max(g.zF[1], min(g.zF[end], phyt[3]))
end

"""
    agent_diffusionX(phyt,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionX(phyt,g,κh)
    phyt[1] += rand(Uniform(-1.0,1.0)) * κh
    phyt[1] = periodic_domain(g.xF, phyt[1])
end

"""
    agent_diffusionX(phyt,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionY(phyt,g,κh)
    phyt[2] += rand(Uniform(-1.0,1.0)) * κh
    phyt[2] = periodic_domain(g.yF, phyt[2])
end
"""
    agent_diffusionZ(phyt,g,κv)
Using a random walk algorithm for vertical diffusion
"""
function agent_diffusionZ(phyt,g,κv)
    phyt[3] += rand(Uniform(-1.0,1.0)) * κv
    phyt[3] = max(g.zF[1], min(g.zF[end], phyt[3]))
end
