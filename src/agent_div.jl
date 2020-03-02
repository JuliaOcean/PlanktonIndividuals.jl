#########################################################################
# advection  of agents (use Interpolations.jl to interpolate velocity fields)#
#########################################################################
"""
    generate_vel_itp(grid, vel)
Use Interpolations.jl to generate interpolation objects
'Ogrid' is the grid information of Oceananigans.jl
"""
function generate_vel_itp(Ogrid, vel)
    # deal with halo points
    xF = collect(Ogrid.xF)
    pushfirst!(xF,xF[1]-(xF[2]-xF[1]))
    yF = collect(Ogrid.yF)
    pushfirst!(yF,yF[1]-(yF[2]-yF[1]))
    zF = collect(Ogrid.zF)
    pushfirst!(zF,zF[1]-(zF[2]-zF[1]))
    pushfirst!(zF,zF[1]-(zF[2]-zF[1]))

    xC = collect(Ogrid.xC)
    pushfirst!(xC,xF[2]-(xC[1]-xF[2]))
    push!(xC,xF[end]+(xF[end]-xC[end]))
    yC = collect(Ogrid.yC)
    pushfirst!(yC,yF[2]-(yC[1]-yF[2]))
    push!(yC,yF[end]+(yF[end]-yC[end]))
    zC = collect(Ogrid.zC)
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
    uvel, vvel, wvel = get_vels(phyt.x, phyt.y, phyt.z, vel_itp)
    phyt.x = phyt.x + uvel*ΔT
    phyt.y = phyt.y + vvel*ΔT
    phyt.z = max(g.zF[1],min(g.zF[end],phyt.z + wvel*ΔT))
    # periodic domain
    phyt.x = periodic_domain(g.xF, phyt.x)
    phyt.y = periodic_domain(g.yF, phyt.y)
end
"""
    agent_advectionRK4(phyts_a, vel_itps, g, ΔT::Int64)
'vel_itps' is an array of tuples containing interpolations of u, v, w velocites of current time step
"""
function agent_advectionRK4(phyt, vel_itps, g, ΔT::Int64)
    u1,v1,w1 = get_vels(phyt.x, phyt.y, phyt.z, vel_itps[1]) # velocites at t
    gx1 = periodic_domain(g.xF, phyt.x + u1*0.5*ΔT)
    gy1 = periodic_domain(g.yF, phyt.y + v1*0.5*ΔT)
    gz1 = phyt.z + w1*0.5*ΔT
    gz1 = max(g.zF[1],min(g.zF[end],gz1))
    u2,v2,w2 = get_vels(gx1, gy1, gz1, vel_itps[2]) # velocites at t+0.5ΔT
    gx2 = periodic_domain(g.xF, phyt.x + u2*0.5*ΔT)
    gy2 = periodic_domain(g.yF, phyt.y + v2*0.5*ΔT)
    gz2 = phyt.z + w2*0.5*ΔT
    gz2 = max(g.zF[1],min(g.zF[end],gz2))
    u3,v3,w3 = get_vels(gx2, gy2, gz2, vel_itps[2]) # velocites at t+0.5ΔT
    gx3 = periodic_domain(g.xF, phyt.x + u3*0.5*ΔT)
    gy3 = periodic_domain(g.yF, phyt.y + v3*0.5*ΔT)
    gz3 = phyt.z + w3*0.5*ΔT
    gz3 = max(g.zF[1],min(g.zF[end],gz3))
    u4,v4,w4 = get_vels(gx3, gy3, gz3, vel_itps[3]) # velocites at t+ΔT
    dx = (u1 + 2*u2 + 2*u3 + u4) / 6 * ΔT
    dy = (v1 + 2*v2 + 2*v3 + v4) / 6 * ΔT
    dz = (w1 + 2*w2 + 2*w3 + w4) / 6 * ΔT
    phyt.x = periodic_domain(g.xF, phyt.x + dx)
    phyt.y = periodic_domain(g.yF, phyt.y + dy)
    phyt.z = phyt.z + dz
    phyt.z = max(g.zF[1], min(g.zF[end], phyt.z))
end

"""
    agent_diffusionX(phyt,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionX(phyt,g,κh)
    phyt.x += rand(Uniform(-1.0,1.0)) * κh
    phyt.x = periodic_domain(g.xF, phyt.x)
end

"""
    agent_diffusionX(phyt,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionY(phyt,g,κh)
    phyt.y += rand(Uniform(-1.0,1.0)) * κh
    phyt.y = periodic_domain(g.yF, phyt.y)
end
"""
    agent_diffusionZ(phyt,g,κv)
Using a random walk algorithm for vertical diffusion
"""
function agent_diffusionZ(phyt,g,κv)
    phyt.z += rand(Uniform(-1.0,1.0)) * κv
    phyt.z = max(g.zF[1], min(g.zF[end], phyt.z))
end

# """
#     simple_itpl(x, y, z, a)
# simple interpolation: interpolate according to c grid (velocity on faces), for 1d only
# 'x', 'y', 'z' are grid indices, 'a' is the w velocity field need to interpolate
# """
# function simple_itpl(x, y, z, a)
#     x₀, y₀, z₀ = trunc(int,x), trunc(int,y), trunc(int,z)
#     zᵈ = z - z₀
#     w₋ = a[x₀, y₀, z₀]
#     w₊ = a[x₀, y₀, z₀+1]
#     wvel = w₋ * (1 - zᵈ) + w₊ * zᵈ
#     return wvel
# end

# """
#     bilinear_itlp(x, y, z, a)
# bilinear interpolation of horizontal velocities
# 'x', 'y', 'z' are grid indices, 'a' is the velocity field need to interpolate, e.g. u, v
# """
# function bilinear_itpl(x, y, z, a)
#     x₀, y₀, z₀ = trunc(int,x), trunc(int,y), trunc(int,z)
#     xᵈ = x - x₀
#     yᵈ = y - y₀
#     vel_00 = a[x₀, y₀, z₀]
#     vel_10 = a[x₀+1, y₀, z₀]
#     vel_01 = a[x₀, y₀+1, z₀]
#     vel_11 = a[x₀+1, y₀+1, z₀]
#     vel_0 = vel_00 * (1 - yᵈ) + vel_10 * yᵈ
#     vel_1 = vel_01 * (1 - yᵈ) + vel_11 * yᵈ
#     vel = vel_0 * (1 - xᵈ) + vel_1 * xᵈ
#     return vel
# end

# """
#     double_grid_2D(velᵇ,grid)
# Compute double grid of each time step
# 'velᵇ' is the velocitiy fields on big grids of current time step
# 'velᵇ' is required to have extra cols & rows at the head and bottom of each dim
# """
# function double_grid_2D(velᵇ)
#     Nx, Ny, Nz = size(velᵇ.u)
#     ny2h = Ny*2-1; nx2h = Nx*2-1;
#     vel_sh = (nx2h, ny2h, Nz);
#     u2h = zeros(vel_sh);
#     v2h = zeros(vel_sh);
#     # Compute values for new grids
#     uX = 0.5*(velᵇ.u[1:Nx-1,:,:] + velᵇ.u[2:Nx,:,:]);
#     vY = 0.5*(velᵇ.v[:,1:Ny-1,:] + velᵇ.v[:,2:Ny,:]);
#     u2h[2:2:nx2h,1:2:ny2h,:] = uX;
#     u2h[1:2:nx2h,1:2:ny2h,:] = velᵇ.u;
#     v2h[1:2:nx2h,2:2:ny2h,:] = vY;
#     v2h[1:2:nx2h,1:2:ny2h,:] = velᵇ.v;
#     uY = 0.5*(u2h[:,1:2:ny2h-1,:] + u2h[:,3:2:ny2h,:]);
#     vX = 0.5*(v2h[1:2:nx2h-1,:,:] + v2h[3:2:nx2h,:,:]);
#     u2h[:,2:2:ny2h,:] = uY;
#     v2h[2:2:nx2h,:,:] = vX;
#     # Delete the boundaries
#     uvel = u2h[3:end,2:end-1,:];
#     vvel = v2h[2:end-1,3:end,:];
#     # Deal with vertical velocities
#     if Nz > 1
#         wvel = zeros(vel_sh)
#         wvel[1:2:nx2h-1,1:2:ny2h-1,:] = velᵇ.w[1:Nx-1,1:Ny-1,:];
#         wvel[1:2:nx2h-1,2:2:ny2h,:] = velᵇ.w[1:Nx-1,1:Ny-1,:];
#         wvel[2:2:nx2h,:,:] = wvel[1:2:nx2h-1,:,:];
#     else
#         wvel = zeros(vel_sh)
#     end
#     velᵈ = velocity(uvel, vvel, wvel)
#     return velᵈ
# end

# """
#     double_grid_3D(velᵇ)
# Compute double grid of each time step from Oceananigans grid (with halo point)
# 'velᵇ' is the velocitiy fields on big grids of current time step
# """
# function double_grid_3D(velᵇ)
#     Nx, Ny, Nz = size(velᵇ.u)
#     ny2h = Ny*2-1; nx2h = Nx*2-1; nz2h = Nz*2-1;
#     vel_sh = (nx2h, ny2h, nz2h);
#     u2h = zeros(vel_sh);
#     v2h = zeros(vel_sh);
#     w2h = zeros(vel_sh);
#     # Compute values for new grids
#     uX = 0.5*(velᵇ.u[1:Nx-1,:,:] + velᵇ.u[2:Nx,:,:]);
#     vY = 0.5*(velᵇ.v[:,1:Ny-1,:] + velᵇ.v[:,2:Ny,:]);
#     wZ = 0.5*(velᵇ.w[:,:,1:Nz-1] + velᵇ.w[:,:,2:Nz]);
#     u2h[2:2:nx2h,1:2:ny2h,1:2:nz2h] = uX;
#     u2h[1:2:nx2h,1:2:ny2h,1:2:nz2h] = velᵇ.u;
#     v2h[1:2:nx2h,2:2:ny2h,1:2:nz2h] = vY;
#     v2h[1:2:nx2h,1:2:ny2h,1:2:nz2h] = velᵇ.v;
#     w2h[1:2:nx2h,1:2:ny2h,2:2:nz2h] = wZ;
#     w2h[1:2:nx2h,1:2:ny2h,1:2:nz2h] = velᵇ.w;
#     uY = 0.5*(u2h[:,1:2:ny2h-1,:] + u2h[:,3:2:ny2h,:]);
#     vX = 0.5*(v2h[1:2:nx2h-1,:,:] + v2h[3:2:nx2h,:,:]);
#     wX = 0.5*(w2h[1:2:nx2h-1,:,:] + w2h[3:2:nx2h,:,:]);
#     u2h[:,2:2:ny2h,:] = uY;
#     v2h[2:2:nx2h,:,:] = vX;
#     w2h[2:2:nx2h,:,:] = wX;
#     uZ = 0.5*(u2h[:,:,1:2:nz2h-1] + u2h[:,:,3:2:nz2h]);
#     vZ = 0.5*(v2h[:,:,1:2:nz2h-1] + v2h[:,:,3:2:nz2h]);
#     wY = 0.5*(w2h[:,1:2:ny2h-1,:] + w2h[:,3:2:ny2h,:]);
#     u2h[:,:,2:2:nz2h] = uZ;
#     v2h[:,:,2:2:nz2h] = vZ;
#     w2h[:,2:2:nz2h,:] = wY;
#     # Delete the boundaries
#     uvel = u2h[3:end,2:end-1,2:end-1];
#     vvel = v2h[2:end-1,3:end,2:end-1];
#     wvel = w2h[2:end-1,2:end-1,3:end];
#     velᵈ = velocity(uvel, vvel, wvel);
#     return velᵈ
# end
# """
#     trilinear_itlp(x, y, z, a)
# Trilinear interpolation of velocities
# 'x', 'y', 'z' are doubled grid indices(whether 2D or 3D),
# 'a' is the velocity field need to interpolate, e.g. u, v, w
# """
# function trilinear_itpl(x, y, z, a)
#     x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
#     xᵈ = x - x₀
#     yᵈ = y - y₀
#     zᵈ = z - z₀
#     vel_000 = a[x₀, y₀, z₀]
#     vel_100 = a[x₀+1, y₀, z₀]
#     vel_001 = a[x₀, y₀, z₀+1]
#     vel_010 = a[x₀, y₀+1, z₀]
#     vel_110 = a[x₀+1, y₀+1, z₀]
#     vel_011 = a[x₀, y₀+1, z₀+1]
#     vel_101 = a[x₀+1, y₀, z₀+1]
#     vel_111 = a[x₀+1, y₀+1, z₀+1]
#     vel_00 = vel_000 * (1 - yᵈ) + vel_100 * yᵈ
#     vel_01 = vel_001 * (1 - yᵈ) + vel_101 * yᵈ
#     vel_10 = vel_010 * (1 - yᵈ) + vel_110 * yᵈ
#     vel_11 = vel_011 * (1 - yᵈ) + vel_111 * yᵈ
#     vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
#     vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
#     vel = vel_0 * (1-zᵈ) + vel_1 * zᵈ
#     return vel
# end
# """
#     get_vels(x, y, z, g, velᵈ)
# Read velocities at (x,y,z) from velocity fields
# 'x', 'y', 'z' are original grid indices
# the velocity field passed to the function can be boubled grids or original grids
# """
# function get_vels(x, y, z, g, vels, grid_type::String)
#     if grid_type == "3D"
#         uvel = trilinear_itpl(2*x-1, 2*y-1, 2*z-1, vels.u)
#         vvel = trilinear_itpl(2*x-1, 2*y-1, 2*z-1, vels.v)
#         wvel = trilinear_itpl(2*x-1, 2*y-1, 2*z-1, vels.w)
#     elseif grid_type == "2D"
#         if g.Nz >1
#             g.Nx == 1 ? uvel = 0.0 : uvel = trilinear_itpl(2*x-1, 2*y-1, z, vels.u) # unit: m/s, trilinear interpolation
#             g.Ny == 1 ? vvel = 0.0 : vvel = trilinear_itpl(2*x-1, 2*y-1, z, vels.v) # unit: m/s, trilinear interpolation
#         else
#             g.Nx == 1 ? uvel = 0.0 : uvel = bilinear_itpl(2*x-1, 2*y-1, z, vels.u) # unit: m/s, bilinear interpolation
#             g.Ny == 1 ? vvel = 0.0 : vvel = bilinear_itpl(2*x-1, 2*y-1, z, vels.v) # unit: m/s, bilinear interpolation
#         end
#         wvel = trilinear_itpl(2*x-1, 2*y-1, z, vels.w) # unit: m/s, trilinear interpolation
#     elseif grid_type == "1D"
#         uvel = 0.0; vvel = 0.0
#         wvel = simple_itpl(x,y,z, vels.w) # unit: m/s, simple interpolation
#     else
#         return "'grid_type' should be '3D', '2D' or '1D'"
#     end
#     return uvel, vvel, wvel
# end
# trilinear interpolation: based on C grid(local double gird,velocity on corners)
#function trilinear_local_itpl(x,y,z,velocity)
#    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
#    vel = zeros(3,3,2,3)
#    vel[1,1,1,1] = 0.5 * (velocity.u[x₀,y₀,z₀] + velocity.u[x₀,y₀-1,z₀])
#    vel[2,1,1,1] = velocity.u[x₀,y₀,z₀]
#    vel[3,1,1,1] = 0.5 * (velocity.u[x₀,y₀,z₀] + velocity.u[x₀,y₀+1,z₀])
#    vel[1,3,1,1] = 0.5 * (velocity.u[x₀+1,y₀,z₀] + velocity.u[x₀+1,y₀-1,z₀])
#    vel[2,3,1,1] = velocity.u[x₀+1,y₀,z₀]
#    vel[3,3,1,1] = 0.5 * (velocity.u[x₀+1,y₀,z₀] + velocity.u[x₀+1,y₀+1,z₀])
#    vel[1,2,1,1] = 0.5 * (vel[1,1,1,1] + vel[1,3,1,1])
#    vel[2,2,1,1] = 0.5 * (vel[2,1,1,1] + vel[2,3,1,1])
#    vel[3,2,1,1] = 0.5 * (vel[3,1,1,1] + vel[3,3,1,1])
#
#    vel[1,1,2,1] = 0.5 * (velocity.u[x₀,y₀,z₀+1] + velocity.u[x₀,y₀-1,z₀+1])
#    vel[2,1,2,1] = velocity.u[x₀,y₀,z₀+1]
#    vel[3,1,2,1] = 0.5 * (velocity.u[x₀,y₀,z₀+1] + velocity.u[x₀,y₀+1,z₀+1])
#    vel[1,3,2,1] = 0.5 * (velocity.u[x₀+1,y₀,z₀+1] + velocity.u[x₀+1,y₀-1,z₀+1])
#    vel[2,3,2,1] = velocity.u[x₀+1,y₀,z₀+1]
#    vel[3,3,2,1] = 0.5 * (velocity.u[x₀+1,y₀,z₀+1] + velocity.u[x₀+1,y₀+1,z₀+1])
#    vel[1,2,2,1] = 0.5 * (vel[1,1,2,1] + vel[1,3,2,1])
#    vel[2,2,2,1] = 0.5 * (vel[2,1,2,1] + vel[2,3,2,1])
#    vel[3,2,2,1] = 0.5 * (vel[3,1,2,1] + vel[3,3,2,1])
#
#    vel[1,1,1,2] = 0.5 * (velocity.v[x₀-1,y₀,z₀] + velocity.v[x₀,y₀,z₀])
#    vel[1,2,1,2] = velocity.v[x₀,y₀,z₀]
#    vel[1,3,1,2] = 0.5 * (velocity.v[x₀,y₀,z₀] + velocity.v[x₀+1,y₀,z₀])
#    vel[3,1,1,2] = 0.5 * (velocity.v[x₀-1,y₀+1,z₀] + velocity.v[x₀,y₀+1,z₀])
#    vel[3,2,1,2] = velocity.v[x₀,y₀+1,z₀]
#    vel[3,3,1,2] = 0.5 * (velocity.v[x₀,y₀+1,z₀] + velocity.v[x₀+1,y₀+1,z₀])
#    vel[2,1,1,2] = 0.5 * (vel[1,1,1,2] + vel[3,1,1,2])
#    vel[2,2,1,2] = 0.5 * (vel[1,2,1,2] + vel[3,2,1,2])
#    vel[2,3,1,2] = 0.5 * (vel[1,3,1,2] + vel[3,3,1,2])
#
#    vel[1,1,2,2] = 0.5 * (velocity.v[x₀-1,y₀,z₀+1] + velocity.v[x₀,y₀,z₀+1])
#    vel[1,2,2,2] = velocity.v[x₀,y₀,z₀+1]
#    vel[1,3,2,2] = 0.5 * (velocity.v[x₀,y₀,z₀+1] + velocity.v[x₀+1,y₀,z₀+1])
#    vel[3,1,2,2] = 0.5 * (velocity.v[x₀-1,y₀+1,z₀+1] + velocity.v[x₀,y₀+1,z₀+1])
#    vel[3,2,2,2] = velocity.v[x₀,y₀+1,z₀+1]
#    vel[3,3,2,2] = 0.5 * (velocity.v[x₀,y₀+1,z₀+1] + velocity.v[x₀+1,y₀+1,z₀+1])
#    vel[2,1,2,2] = 0.5 * (vel[1,1,2,2] + vel[3,1,2,2])
#    vel[2,2,2,2] = 0.5 * (vel[1,2,2,2] + vel[3,2,2,2])
#    vel[2,3,2,2] = 0.5 * (vel[1,3,2,2] + vel[3,3,2,2])
#
#    vel[1,1,1,3] = velocity.w[x₀-1,y₀-1,z₀]
#    vel[1,2,1,3] = vel[1,3,1,3] = velocity.w[x₀,y₀-1,z₀]
#    vel[2,1,1,3] = vel[3,1,1,3] = velocity.w[x₀-1,y₀,z₀]
#    vel[2,2,1,3] = vel[3,2,1,3] = vel[2,3,1,3] = vel[3,3,1,3] = velocity.w[x₀,y₀,z₀]
#
#    vel[1,1,2,3] = velocity.w[x₀-1,y₀-1,z₀+1]
#    vel[1,2,2,3] = vel[1,3,2,3] = velocity.w[x₀,y₀-1,z₀+1]
#    vel[2,1,2,3] = vel[3,1,2,3] = velocity.w[x₀-1,y₀,z₀+1]
#    vel[2,2,2,3] = vel[3,2,2,3] = vel[2,3,2,3] = vel[3,3,2,3] = velocity.w[x₀,y₀,z₀+1]
#
#    xᵈ = x - x₀
#    yᵈ = y - y₀
#    zᵈ = z - z₀
#    if xᵈ ≤ 0.5 && yᵈ ≤ 0.5
#        vel_00 = vel[1,1,1,:] .* (1 - yᵈ) .+ vel[2,1,1,:] .* yᵈ
#        vel_01 = vel[1,1,2,:] .* (1 - yᵈ) .+ vel[2,1,2,:] .* yᵈ
#        vel_10 = vel[1,2,1,:] .* (1 - yᵈ) .+ vel[2,2,1,:] .* yᵈ
#        vel_11 = vel[1,2,2,:] .* (1 - yᵈ) .+ vel[2,2,2,:] .* yᵈ
#        vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
#        vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
#        vel₀ = vel_0 * (1-zᵈ) + vel_1 * zᵈ
#    elseif xᵈ ≤ 0.5 && yᵈ > 0.5
#        vel_00 = vel[2,1,1,:] .* (1 - yᵈ) .+ vel[3,1,1,:] .* yᵈ
#        vel_01 = vel[2,1,2,:] .* (1 - yᵈ) .+ vel[3,1,2,:] .* yᵈ
#        vel_10 = vel[2,2,1,:] .* (1 - yᵈ) .+ vel[3,2,1,:] .* yᵈ
#        vel_11 = vel[2,2,2,:] .* (1 - yᵈ) .+ vel[3,2,2,:] .* yᵈ
#        vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
#        vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
#        vel₀ = vel_0 * (1-zᵈ) + vel_1 * zᵈ
#    elseif xᵈ > 0.5 && yᵈ ≤ 0.5
#        vel_00 = vel[1,2,1,:] .* (1 - yᵈ) .+ vel[2,2,1,:] .* yᵈ
#        vel_01 = vel[1,2,2,:] .* (1 - yᵈ) .+ vel[2,2,2,:] .* yᵈ
#        vel_10 = vel[1,3,1,:] .* (1 - yᵈ) .+ vel[2,3,1,:] .* yᵈ
#        vel_11 = vel[1,3,2,:] .* (1 - yᵈ) .+ vel[2,3,2,:] .* yᵈ
#        vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
#        vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
#        vel₀ = vel_0 * (1-zᵈ) + vel_1 * zᵈ
#    else # xᵈ > 0.5 && yᵈ > 0.5
#        vel_00 = vel[2,2,1,:] .* (1 - yᵈ) .+ vel[3,2,1,:] .* yᵈ
#        vel_01 = vel[2,2,2,:] .* (1 - yᵈ) .+ vel[3,2,2,:] .* yᵈ
#        vel_10 = vel[2,3,1,:] .* (1 - yᵈ) .+ vel[3,3,1,:] .* yᵈ
#        vel_11 = vel[2,3,2,:] .* (1 - yᵈ) .+ vel[3,3,2,:] .* yᵈ
#        vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
#        vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
#        vel₀ = vel_0 * (1-zᵈ) + vel_1 * zᵈ
#    end
#    return vel₀
#end
