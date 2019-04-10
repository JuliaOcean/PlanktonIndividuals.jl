# Increment and decrement integer a with periodic wrapping.
@inline incmod1(a, n) = ifelse(a==n+1, 1, a + 1)
@inline decmod1(a, n) = ifelse(a==1, n, a - 1)
# i represents longitude, j represents latitude
# in Julia array, count rows first then columns, so it should be a[j, i, k] to represent a grid.
@inline avgx_c2f(g::grids, f, i, j, k) = @inbounds 0.5 * (f[j, i, k] + f[j, decmod1(i, g.Nx), k])
@inline avgy_c2f(g::grids, f, i, j, k) = @inbounds 0.5 * (f[j, i, k] + f[decmod1(j, g.Ny), i, k])
@inline function avgz_c2f(g::grids, f, i, j, k)
    if k == 1
        @inbounds return f[j, i, k]
    else
        @inbounds return  0.5 * (f[j, i, k] + f[j, i, k-1])
    end
end
@inline δx_c2f(g::grids, f, i, j, k) = @inbounds f[j, i, k] - f[j, decmod1(i, g.Nx), k]
@inline δy_c2f(g::grids, f, i, j, k) = @inbounds f[j, i, k] - f[decmod1(j, g.Ny), i, k]
@inline function δz_c2f(g::grids, f, i, j, k)
    if k == 1
        return 0
    else
        @inbounds return f[j, i, k-1] - f[j, i, k]
    end
end
@inline function δx_f2c_ab̄ˣ(g::grids, a, b, i, j, k)
    @inbounds (g.Δy[j] * g.Δz[k] * a[j, incmod1(i, g.Nx), k] * avgx_c2f(g, b, incmod1(i, g.Nx), j, k) - 
               g.Δy[j] * g.Δz[k] * a[j, i, k] * avgx_c2f(g, b, i, j, k))
end

@inline function δy_f2c_ab̄ʸ(g::grids, a, b, i, j, k)
    @inbounds (g.Δx[i] * g.Δz[k] * a[incmod1(j, g.Ny), i, k] * avgy_c2f(g, b, i, incmod1(j, g.Ny), k) - 
               g.Δx[i] * g.Δz[k] * a[j, i, k] * avgy_c2f(g, b, i, j, k))
end

@inline function δz_f2c_ab̄ᶻ(g::grids, a, b, i, j, k)
    if k == g.Nz
        @inbounds return g.Δx[i] * g.Δy[j] * a[j, i, k] * avgz_c2f(g, b, i, j, k)
    else
        @inbounds return (g.Δx[i] * g.Δy[j] * a[j, i, k] * avgz_c2f(g, b, i, j, k) - 
                          g.Δx[i] * g.Δy[j] * a[j, i, k+1] * avgz_c2f(g, b, i, j, k+1))
    end
end

@inline function div_flux(g::grids, u, v, w, Q, i, j, k)
    @inbounds ΔV = g.Δx[i] * g.Δy[j] * g.Δz[k]
    if k == 1
        @inbounds return (δx_f2c_ab̄ˣ(g, u, Q, i, j, k) + δy_f2c_ab̄ʸ(g, v, Q, i, j, k) - w[j, i, 2] * avgz_c2f(g, Q, i, j, 2)) / ΔV 
    else
        return (δx_f2c_ab̄ˣ(g, u, Q, i, j, k) + δy_f2c_ab̄ʸ(g, v, Q, i, j, k) + δz_f2c_ab̄ᶻ(g, w, Q, i, j, k)) / ΔV
    end
end

@inline function δx²_c2f2c(g::grids, f, i, j, k)
    @inbounds (κh * g.Δy[j] * g.Δz[k] * δx_c2f(g, f, incmod1(i, g.Nx), j, k) - 
               κh * g.Δy[j] * g.Δz[k] * δx_c2f(g, f, i, j, k))
end
@inline function δy²_c2f2c(g::grids, f, i, j, k)
    @inbounds (κh * g.Δx[i] * g.Δz[k] * δy_c2f(g, f, i, incmod1(j, g.Ny), k) - 
               κh * g.Δx[i] * g.Δz[k] * δy_c2f(g, f, i, j, k))
end
@inline function δz²_c2f2c(g::grids, f, i, j, k)
    if k == g.Nz
        return κv * g.Δx[i] * g.Δy[j] * δz_c2f(g, f, i, j, k)
    else
        return (κv * g.Δx[i] * g.Δy[j] * δz_c2f(g, f, i, j, k) - 
                κv * g.Δx[i] * g.Δy[j] * δz_c2f(g, f, i, j, k+1))
    end
end
@inline function κ∇²(g::grids, Q, κh, κv, i, j, k)
    @inbounds ΔV = g.Δx[i] * g.Δy[j] * g.Δz[k]
    return (δx²_c2f2c(g, Q, i, j, k) + δy²_c2f2c(g, Q, i, j, k) + δz²_c2f2c(g, Q, i, j, k)) / ΔV
end
