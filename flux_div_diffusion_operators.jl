# Increment and decrement integer a with periodic wrapping.
@inline incmod1(a, n) = ifelse(a==n+1, 1, a + 1)
@inline decmod1(a, n) = ifelse(a==1, n, a - 1)
@inline avgx_c2f(g::grids, f, i, j, k) = @inbounds 0.5 * (f[i, j, k] + f[decmod1(i, g.Nx), j, k])
@inline avgy_c2f(g::grids, f, i, j, k) = @inbounds 0.5 * (f[i, j, k] + f[i, decmod1(j, g.Ny), k])
@inline function avgz_c2f(g::grids, f, i, j, k)
    if k == 1
        @inbounds return f[i, j, k]
    else
        @inbounds return  0.5 * (f[i, j, k] + f[i, j, k-1])
    end
end
@inline δx_c2f(g::grids, f, i, j, k) = @inbounds f[i, j, k] - f[decmod1(i, g.Nx), j, k]
@inline δy_c2f(g::grids, f, i, j, k) = @inbounds f[i, j, k] - f[i, decmod1(j, g.Ny), k]
@inline function δz_c2f(g::grids, f, i, j, k)
    if k == 1
        return 0
    else
        @inbounds return f[i, j, k-1] - f[i, j, k]
    end
end
@inline function δx_f2c_ab̄ˣ(g::grids, a, b, i, j, k)
    @inbounds (g.Δy[incmod1(i, g.Nx), j, k] * g.Δz[incmod1(i, g.Nx), j, k] * a[incmod1(i, g.Nx), j, k] * avgx_c2f(g, b, incmod1(i, g.Nx), j, k) - 
               g.Δy[i, j, k] * g.Δz[i, j, k] * a[i, j, k] * avgx_c2f(g, b, i, j, k))
end

@inline function δy_f2c_ab̄ʸ(g::grids, a, b, i, j, k)
    @inbounds (g.Δx[i, incmod1(j, g.Ny), k] * g.Δz[i, incmod1(j, g.Ny), k] * a[i, incmod1(j, g.Ny), k] * avgy_c2f(g, b, i, incmod1(j, g.Ny), k) - 
               g.Δx[i, j, k] * g.Δz[i, j, k] * a[i, j, k] * avgy_c2f(g, b, i, j, k))
end

@inline function δz_f2c_ab̄ᶻ(g::grids, a, b, i, j, k)
    if k == g.Nz
        @inbounds return g.Δx[i, j, k] * g.Δy[i, j, k] * a[i, j, k] * avgz_c2f(g, b, i, j, k)
    else
        @inbounds return (g.Δx[i, j, k] * g.Δy[i, j, k] * a[i, j, k] * avgz_c2f(g, b, i, j, k) - 
                          g.Δx[i, j, k+1] * g.Δy[i, j, k+1] * a[i, j, k+1] * avgz_c2f(g, b, i, j, k+1))
    end
end

@inline function div_flux(g::grids, u, v, w, Q, i, j, k)
    @inbounds ΔV = g.Δx[i, j, k] * g.Δy[i, j, k] * g.Δz[i, j, k]
    if k == 1
        @inbounds return (δx_f2c_ab̄ˣ(g, u, Q, i, j, k) + δy_f2c_ab̄ʸ(g, v, Q, i, j, k) - w[i, j, 2] * avgz_c2f(g, Q, i, j, 2)) / ΔV 
    else
        return (δx_f2c_ab̄ˣ(g, u, Q, i, j, k) + δy_f2c_ab̄ʸ(g, v, Q, i, j, k) + δz_f2c_ab̄ᶻ(g, w, Q, i, j, k)) / ΔV
    end
end

@inline function δx²_c2f2c(g::grids, f, i, j, k)
    @inbounds (κh * g.Δy[incmod1(i, g.Nx), j, k] * g.Δz[incmod1(i, g.Nx), j, k] * δx_c2f(g, f, incmod1(i, g.Nx), j, k) - 
               κh * g.Δy[i, j, k] * g.Δz[i, j, k] * δx_c2f(g, f, i, j, k))
end
@inline function δy²_c2f2c(g::grids, f, i, j, k)
    @inbounds (κh * g.Δx[i, incmod1(j, g.Ny), k] * g.Δz[i, incmod1(j, g.Ny), k] * δy_c2f(g, f, i, incmod1(j, g.Ny), k) - 
               κh * g.Δx[i, j, k] * g.Δz[i, j, k] * δy_c2f(g, f, i, j, k))
end
@inline function δz²_c2f2c(g::grids, f, i, j, k)
    if k == g.Nz
        return κv * g.Δx[i, j, k] * g.Δy[i, j, k] * δz_c2f(g, f, i, j, k)
    else
        return (κv * g.Δx[i, j, k] * g.Δy[i, j, k] * δz_c2f(g, f, i, j, k) - 
                κv * g.Δx[i, j, k+1] * g.Δy[i, j, k+1] * δz_c2f(g, f, i, j, k+1))
    end
end
@inline function κ∇²(g::grids, Q, κh, κv, i, j, k)
    @inbounds ΔV = g.Δx[i, j, k] * g.Δy[i, j, k] * g.Δz[i, j, k]
    return (δx²_c2f2c(g, Q, i, j, k) + δy²_c2f2c(g, Q, i, j, k) + δz²_c2f2c(g, Q, i, j, k)) / ΔV
end
