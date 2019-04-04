# Increment and decrement integer a with periodic wrapping.
@inline incmod1(a, n) = ifelse(a==n, 1, a + 1)
@inline decmod1(a, n) = ifelse(a==1, n, a - 1)
@inline avgx_c2f(g::grid, f, i, j, k) = @inbounds 0.5 * (f[i, j, k] + f[decmod1(i, g.Nx), j, k])
@inline avgy_c2f(g::grid, f, i, j, k) = @inbounds 0.5 * (f[i, j, k] + f[i, decmod1(j, g.Ny), k])
@inline function avgz_c2f(g::grid, f, i, j, k)
    if k == 1
        @inbounds return f[i, j, k]
    else
        @inbounds return  0.5 * (f[i, j, k] + f[i, j, k-1])
    end
end

@inline function δx_f2c_ab̄ˣ(g::grid, a, b, i, j, k)
    @inbounds (a[incmod1(i, g.Nx), j, k] * avgx_c2f(g, b, incmod1(i, g.Nx), j, k) - a[i,j, k] * avgx_c2f(g, b, i,j, k))
end

@inline function δy_f2c_ab̄ʸ(g::grid, a, b, i, j, k)
    @inbounds (a[i, incmod1(j, g.Ny), k] * avgy_c2f(g, b, i, incmod1(j, g.Ny), k) - a[i,j, k] * avgy_c2f(g, b, i, j,k))
end

@inline function δz_f2c_ab̄ᶻ(g::grid, a, b, i, j, k)
    if k == g.Nz
        @inbounds return a[i, j, k] * avgz_c2f(g, b, i, j, k)
    else
        @inbounds return (a[i, j,   k] * avgz_c2f(g, b, i, j,   k) - a[i, j, k+1] * avgz_c2f(g, b, i, j, k+1))
    end
end

@inline function div_flux(g::grid, u, v, w, Q, i, j, k)
    if k == 1
        @inbounds return (δx_f2c_ab̄ˣ(g, u, Q, i, j, k) / g.Δx) + (δy_f2c_ab̄ʸ(g, v, Q, i, j, k) / g.Δy) - ((w[i, j, 2] * avgz_c2f(g, Q, i, j, 2)) / g.Δz)
    else
        return (δx_f2c_ab̄ˣ(g, u, Q, i, j, k) / g.Δx) + (δy_f2c_ab̄ʸ(g, v, Q, i, j, k) / g.Δy) + (δz_f2c_ab̄ᶻ(g, w, Q, i, j, k) / g.Δz)
    end
end
