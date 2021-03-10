##### convert the coordinates into fractional indices at grid cell center and face for a individual at (x,y,z) 
##### 0-based indexing and excluding halo regions
@inline get_xc_index(x, grid::Grids) = @inbounds (x - grid.xC[grid.Hx+1]) / grid.Δx 
@inline get_xf_index(x, grid::Grids) = @inbounds (x - grid.xF[grid.Hx+1]) / grid.Δx 

@inline get_yc_index(y, grid::Grids) = @inbounds (y - grid.yC[grid.Hy+1]) / grid.Δy 
@inline get_yf_index(y, grid::Grids) = @inbounds (y - grid.yF[grid.Hy+1]) / grid.Δy 

@inline get_zc_index(z, grid::Grids) = @inbounds (z - grid.zC[grid.Hz+1]) / grid.Δz 
@inline get_zf_index(z, grid::Grids) = @inbounds (z - grid.zF[grid.Hz+1]) / grid.Δz 

##### trilinear interpolation
@inline ψ₀₀₀(xd, yd, zd) = (1 - xd) * (1 - yd) * (1 - zd)
@inline ψ₀₀₁(xd, yd, zd) = (1 - xd) * (1 - yd) *      zd
@inline ψ₀₁₀(xd, yd, zd) = (1 - xd) *      yd  * (1 - zd)
@inline ψ₀₁₁(xd, yd, zd) = (1 - xd) *      yd  *      zd
@inline ψ₁₀₀(xd, yd, zd) =      xd  * (1 - yd) * (1 - zd)
@inline ψ₁₀₁(xd, yd, zd) =      xd  * (1 - yd) *      zd 
@inline ψ₁₁₀(xd, yd, zd) =      xd  *      yd  * (1 - zd)
@inline ψ₁₁₁(xd, yd, zd) =      xd  *      yd  *      zd 

@inline tri_interpolation(vel, xd, yd, zd, xi, yi, zi) = 
        @inbounds (ψ₀₀₀(xd, yd, zd) * vel[xi,   yi,   zi  ] +
                   ψ₀₀₁(xd, yd, zd) * vel[xi,   yi,   zi+1] +
                   ψ₀₁₀(xd, yd, zd) * vel[xi,   yi+1, zi  ] +
                   ψ₀₁₁(xd, yd, zd) * vel[xi,   yi+1, zi+1] +
                   ψ₁₀₀(xd, yd, zd) * vel[xi+1, yi,   zi  ] +
                   ψ₁₀₁(xd, yd, zd) * vel[xi+1, yi,   zi+1] +
                   ψ₁₁₀(xd, yd, zd) * vel[xi+1, yi+1, zi  ] +
                   ψ₁₁₁(xd, yd, zd) * vel[xi+1, yi+1, zi+1])

@inline tri_interpolation_surface(vel, xd, yd, zd, xi, yi, zi) = 
        @inbounds ((ψ₀₀₀(xd, yd, zd) * vel[xi,   yi,   zi] +
                    ψ₀₁₀(xd, yd, zd) * vel[xi,   yi+1, zi] +
                    ψ₁₀₀(xd, yd, zd) * vel[xi+1, yi,   zi] +
                    ψ₁₁₀(xd, yd, zd) * vel[xi+1, yi+1, zi]) / (1 - zd))

@inline function u_itpl(u, x, y, z, ac, g::Grids) 
    xi = get_xf_index(x, g) * ac
    yi = get_yc_index(y, g) * ac
    zi = get_zc_index(z, g) * ac
    xd, xi = mod(xi, 1), unsafe_trunc(Int, xi)
    yd, yi = mod(yi, 1), unsafe_trunc(Int, yi)
    zd, zi = mod(zi, 1), unsafe_trunc(Int, zi)
    ##### shift to 1-based indexing and include halo points
    return tri_interpolation(u, xd, yd, zd, xi+g.Hx+1, yi+g.Hy+1, zi+g.Hz+1)
end

@inline function v_itpl(v, x, y, z, ac, g::Grids) 
    xi = get_xc_index(x, g) * ac
    yi = get_yf_index(y, g) * ac
    zi = get_zc_index(z, g) * ac
    xd, xi = mod(xi, 1), unsafe_trunc(Int, xi)
    yd, yi = mod(yi, 1), unsafe_trunc(Int, yi)
    zd, zi = mod(zi, 1), unsafe_trunc(Int, zi)
    ##### shift to 1-based indexing and include halo points
    return tri_interpolation(v, xd, yd, zd, xi+g.Hx+1, yi+g.Hy+1, zi+g.Hz+1)
end

@inline function w_itpl(w, x, y, z, ac, g::Grids) 
    xi = get_xc_index(x, g) * ac
    yi = get_yc_index(y, g) * ac
    zi = get_zf_index(z, g) * ac
    xd, xi = mod(xi, 1), unsafe_trunc(Int, xi)
    yd, yi = mod(yi, 1), unsafe_trunc(Int, yi)
    zd, zi = mod(zi, 1), unsafe_trunc(Int, zi)
    ##### shift to 1-based indexing and include halo points
    # zi == g.Nz-1 && return tri_interpolation_surface(w, xd, yd, zd, xi+g.Hx+1, yi+g.Hy+1, zi+g.Hz+1)
    return tri_interpolation(w, xd, yd, zd, xi+g.Hx+1, yi+g.Hy+1, zi+g.Hz+1)
end

##### calculate uvw velocities at (x, y, z)
@kernel function vel_interpolate_kernel!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g::Grids)
    i = @index(Global)
    @inbounds uₜ[i] = u_itpl(u, x[i], y[i], z[i], ac[i], g) * ac[i]
    @inbounds vₜ[i] = v_itpl(v, x[i], y[i], z[i], ac[i], g) * ac[i]
    @inbounds wₜ[i] = w_itpl(w, x[i], y[i], z[i], ac[i], g) * ac[i]
end

function vel_interpolate!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g::Grids, arch::Architecture)
    kernel! = vel_interpolate_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g)
    wait(device(arch), event)
    return nothing
end
