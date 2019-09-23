#########################################################################
# advection  of agents (including 3 ways to interpolate velocity fields)#
#########################################################################
# compute double grid of each time step
# velᵇ is the velocitiy fields on big grids of current time step
function double_grid(velᵇ,grid)
    Ny = grid.Ny; Nx = grid.Nx; Nz = grid.Nz;
    u = zeros(Nx*2-1,Ny*2-1,Nz);
    v = zeros(Nx*2-1,Ny*2-1,Nz);
    w = zeros(Nx*2-1,Ny*2-1,Nz);
    velᵈ = velocity(u, v, w);
    ny2h = Ny*2-1; nx2h = Nx*2-1;
    vel_sh = (nx2h, ny2h, Nz);
    u2h = zeros(Float32,vel_sh);
    v2h = zeros(Float32,vel_sh);
    # Compute values for new grids
    uX = 0.5*(velᵇ.u[1:Nx-1,:,:] + velᵇ.u[2:Nx,:,:]);
    vY = 0.5*(velᵇ.v[:,1:Ny-1,:] + velᵇ.v[:,2:Ny,:]);
    u2h[2:2:nx2h,1:2:ny2h,:] = uX;
    u2h[1:2:nx2h,1:2:ny2h,:] = velᵇ.u;
    v2h[1:2:nx2h,2:2:ny2h,:] = vY;
    v2h[1:2:nx2h,1:2:ny2h,:] = velᵇ.v;
    uY = 0.5*(u2h[:,1:2:ny2h-1,:] + u2h[:,3:2:ny2h,:]);
    vX = 0.5*(v2h[1:2:nx2h-1,:,:] + v2h[3:2:nx2h,:,:]);
    u2h[:,2:2:ny2h,:] = uY;
    v2h[2:2:nx2h,:,:] = vX;
    # Delete the boundaries
    velᵈ.u = u2h[2:nx2h,1:ny2h-1,:]; velᵈ.v = v2h[1:nx2h-1,2:ny2h,:];
    # Deal with vertical velocities
    velᵈ.w[1:2:nx2h-1,1:2:ny2h-1,:] = velᵇ.w[1:Nx-1,1:Ny-1,:];
    velᵈ.w[1:2:nx2h-1,2:2:ny2h,:] = velᵇ.w[1:Nx-1,1:Ny-1,:];
    velᵈ.w[2:2:nx2h,:,:] = velᵈ.w[1:2:nx2h-1,:,:];
    return velᵈ
end
function trilinear_itpl(x, y, z, a)
    x = 2*x-2; y=2*y-2;
    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    xᵈ = x - x₀
    yᵈ = y - y₀
    zᵈ = z - z₀
    vel_000 = a[x₀, y₀, z₀]
    vel_100 = a[x₀+1, y₀, z₀]
    vel_001 = a[x₀, y₀, z₀+1]
    vel_010 = a[x₀, y₀+1, z₀]
    vel_110 = a[x₀+1, y₀+1, z₀]
    vel_011 = a[x₀, y₀+1, z₀+1]
    vel_101 = a[x₀+1, y₀, z₀+1]
    vel_111 = a[x₀+1, y₀+1, z₀+1]
    vel_00 = vel_000 * (1 - yᵈ) + vel_100 * yᵈ
    vel_01 = vel_001 * (1 - yᵈ) + vel_101 * yᵈ
    vel_10 = vel_010 * (1 - yᵈ) + vel_110 * yᵈ
    vel_11 = vel_011 * (1 - yᵈ) + vel_111 * yᵈ
    vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
    vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
    vel = vel_0 * (1-zᵈ) + vel_1 * zᵈ
    return vel
end

function agent_move(phyts_a,velᵈ,g,ΔT::Int64)
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        uvel = trilinear_itpl(phyt.x, phyt.y, phyt.z, velᵈ.u) # unit: m/s, trilinear interpolation
        vvel = trilinear_itpl(phyt.x, phyt.y, phyt.z, velᵈ.v) # unit: m/s, trilinear interpolation
        wvel = trilinear_itpl(phyt.x, phyt.y, phyt.z, velᵈ.w) # unit: m/s, trilinear interpolation
#       uvel, vvel, wvel = simple_itpl(phyt.x, phyt.y, phyt.z, vel, t) # unit: m/s, simple interpolation

        xi, yi, zi = trunc(Int,phyt.x), trunc(Int,phyt.y), trunc(Int,phyt.z)
        dx = uvel/g.Lx[xi]*ΔT # unit: grid/h
        dy = vvel/g.Ly[yi]*ΔT # unit: grid/h
        dz = (wvel+k_sink)/g.Lz[zi]*ΔT # vertical movement, plus sinking, unit: grid/h
#       phyt.x = max(1.5,min(g.Nx-0.5,phyt.x - dx*(1+rand()/5)))
#       phyt.y = max(1.5,min(g.Ny-0.5,phyt.y - dy*(1+rand()/5)))
        phyt.x = phyt.x - dx*(1+rand()/3)
        phyt.y = phyt.y - dy*(1+rand()/3)
        phyt.z = max(2.0,min(g.Nz-0.1,phyt.z - dz*(1+rand()/3)))
        # periodic domian
        if phyt.x ≥ g.Nx - 0.5
            phyt.x = phyt.x - g.Nx + 2.0
        end
        if phyt.x ≤ 1.5
            phyt.x = phyt.x + g.Nx - 2.0
        end
        if phyt.y ≥ g.Ny - 0.5
            phyt.y = phyt.y - g.Ny + 2.0
        end
        if phyt.y ≤ 1.5
            phyt.y = phyt.y + g.Ny - 2.0
        end
    end
end
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

# simple interpolation: interpolate according to C grid (velocity on faces), for 1D only
function simple_itpl(x, y, z, vel)
    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    zᵈ = z - z₀
    w₋ = vel.w[x₀, y₀, z₀]
    w₊ = vel.w[x₀, y₀, z₀+1]
    wvel = w₋ * (1 - zᵈ) + w₊ * zᵈ
    return wvel
end
function agent_move_1D(phyts_a,vel,g,ΔT::Int64)
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        wvel = simple_itpl(phyt.x,phyt.y,phyt.z, vel) # unit: m/s, simple interpolation
        zi = trunc(Int,phyt.z)
        dz = wvel/g.Lz[zi]*ΔT # vertical movement, unit: grid/h
        phyt.z = max(1.0,min(g.Nz-0.1,phyt.z - dz*(1+rand()/5)))
    end
end
