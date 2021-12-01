##### convert the coordinates into fractional indices at grid cell center and face for a individual at (x,y,z) 
##### 1-based indexing for grid faces (excluding halo regions)
@inline get_xf_index(x) = x + 1 
@inline get_yf_index(y) = y + 1
@inline get_zf_index(z) = z + 1

##### 1-based indexing for grid centers (excluding halo regions)
@inline get_xc_index(x) = x + 0.5 
@inline get_yc_index(y) = y + 0.5 
@inline get_zc_index(z) = z + 0.5 

##### trilinear interpolation
@inline ψ₀₀₀(xd, yd, zd) = (1 - xd) * (1 - yd) * (1 - zd)
@inline ψ₀₀₁(xd, yd, zd) = (1 - xd) * (1 - yd) *      zd
@inline ψ₀₁₀(xd, yd, zd) = (1 - xd) *      yd  * (1 - zd)
@inline ψ₀₁₁(xd, yd, zd) = (1 - xd) *      yd  *      zd
@inline ψ₁₀₀(xd, yd, zd) =      xd  * (1 - yd) * (1 - zd)
@inline ψ₁₀₁(xd, yd, zd) =      xd  * (1 - yd) *      zd 
@inline ψ₁₁₀(xd, yd, zd) =      xd  *      yd  * (1 - zd)
@inline ψ₁₁₁(xd, yd, zd) =      xd  *      yd  *      zd 

@inline normalized_trilinear_interpolation_u(u, xd, yd, zd, xi, yi, zi, g::AbstractGrid) = 
        @inbounds (ψ₀₀₀(xd, yd, zd) * u[xi,   yi,   zi  ] / ΔxC(xi,   yi,   zi,   g) +
                   ψ₀₀₁(xd, yd, zd) * u[xi,   yi,   zi+1] / ΔxC(xi,   yi,   zi+1, g) +
                   ψ₀₁₀(xd, yd, zd) * u[xi,   yi+1, zi  ] / ΔxC(xi,   yi+1, zi,   g) +
                   ψ₀₁₁(xd, yd, zd) * u[xi,   yi+1, zi+1] / ΔxC(xi,   yi+1, zi+1, g) +
                   ψ₁₀₀(xd, yd, zd) * u[xi+1, yi,   zi  ] / ΔxC(xi+1, yi,   zi,   g) +
                   ψ₁₀₁(xd, yd, zd) * u[xi+1, yi,   zi+1] / ΔxC(xi+1, yi,   zi+1, g) +
                   ψ₁₁₀(xd, yd, zd) * u[xi+1, yi+1, zi  ] / ΔxC(xi+1, yi+1, zi,   g) +
                   ψ₁₁₁(xd, yd, zd) * u[xi+1, yi+1, zi+1] / ΔxC(xi+1, yi+1, zi+1, g))

@inline normalized_trilinear_interpolation_v(v, xd, yd, zd, xi, yi, zi, g::AbstractGrid) = 
        @inbounds (ψ₀₀₀(xd, yd, zd) * v[xi,   yi,   zi  ] / ΔyC(xi,   yi,   zi,   g) +
                   ψ₀₀₁(xd, yd, zd) * v[xi,   yi,   zi+1] / ΔyC(xi,   yi,   zi+1, g) +
                   ψ₀₁₀(xd, yd, zd) * v[xi,   yi+1, zi  ] / ΔyC(xi,   yi+1, zi,   g) +
                   ψ₀₁₁(xd, yd, zd) * v[xi,   yi+1, zi+1] / ΔyC(xi,   yi+1, zi+1, g) +
                   ψ₁₀₀(xd, yd, zd) * v[xi+1, yi,   zi  ] / ΔyC(xi+1, yi,   zi,   g) +
                   ψ₁₀₁(xd, yd, zd) * v[xi+1, yi,   zi+1] / ΔyC(xi+1, yi,   zi+1, g) +
                   ψ₁₁₀(xd, yd, zd) * v[xi+1, yi+1, zi  ] / ΔyC(xi+1, yi+1, zi,   g) +
                   ψ₁₁₁(xd, yd, zd) * v[xi+1, yi+1, zi+1] / ΔyC(xi+1, yi+1, zi+1, g))

@inline normalized_trilinear_interpolation_w(w, xd, yd, zd, xi, yi, zi, g::AbstractGrid) = 
        @inbounds (ψ₀₀₀(xd, yd, zd) * w[xi,   yi,   zi  ] / ΔzC(xi,   yi,   zi,   g) +
                   ψ₀₀₁(xd, yd, zd) * w[xi,   yi,   zi+1] / ΔzC(xi,   yi,   zi+1, g) +
                   ψ₀₁₀(xd, yd, zd) * w[xi,   yi+1, zi  ] / ΔzC(xi,   yi+1, zi,   g) +
                   ψ₀₁₁(xd, yd, zd) * w[xi,   yi+1, zi+1] / ΔzC(xi,   yi+1, zi+1, g) +
                   ψ₁₀₀(xd, yd, zd) * w[xi+1, yi,   zi  ] / ΔzC(xi+1, yi,   zi,   g) +
                   ψ₁₀₁(xd, yd, zd) * w[xi+1, yi,   zi+1] / ΔzC(xi+1, yi,   zi+1, g) +
                   ψ₁₁₀(xd, yd, zd) * w[xi+1, yi+1, zi  ] / ΔzC(xi+1, yi+1, zi,   g) +
                   ψ₁₁₁(xd, yd, zd) * w[xi+1, yi+1, zi+1] / ΔzC(xi+1, yi+1, zi+1, g))

@inline function u_itpl(u, x, y, z, ac, g::AbstractGrid) 
    xi = get_xf_index(x) * ac
    yi = get_yc_index(y) * ac
    zi = get_zc_index(z) * ac
    xd, xi = mod(xi, 1), unsafe_trunc(Int, xi)
    yd, yi = mod(yi, 1), unsafe_trunc(Int, yi)
    zd, zi = mod(zi, 1), unsafe_trunc(Int, zi)
    ##### return fraction of grid spacing
    return normalized_trilinear_interpolation_u(u, xd, yd, zd, xi+g.Hx, yi+g.Hy, zi+g.Hz, g)
end

@inline function v_itpl(v, x, y, z, ac, g::AbstractGrid) 
    xi = get_xc_index(x) * ac
    yi = get_yf_index(y) * ac
    zi = get_zc_index(z) * ac
    xd, xi = mod(xi, 1), unsafe_trunc(Int, xi)
    yd, yi = mod(yi, 1), unsafe_trunc(Int, yi)
    zd, zi = mod(zi, 1), unsafe_trunc(Int, zi)
    ##### return fraction of grid spacing
    return normalized_trilinear_interpolation_v(v, xd, yd, zd, xi+g.Hx, yi+g.Hy, zi+g.Hz, g)
end

@inline function w_itpl(w, x, y, z, ac, g::AbstractGrid) 
    xi = get_xc_index(x) * ac
    yi = get_yc_index(y) * ac
    zi = get_zf_index(z) * ac
    xd, xi = mod(xi, 1), unsafe_trunc(Int, xi)
    yd, yi = mod(yi, 1), unsafe_trunc(Int, yi)
    zd, zi = mod(zi, 1), unsafe_trunc(Int, zi)
    ##### return fraction of grid spacing
    return normalized_trilinear_interpolation_w(w, xd, yd, zd, xi+g.Hx, yi+g.Hy, zi+g.Hz, g)
end