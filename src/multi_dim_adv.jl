##### multi dimensional advection
##### For incompressible model only

##### calculate the divergence of the flux of tracer q only advected by u or v or w
@inline div_flux_x(i, j, k, g::grids, u, q, ΔT) = 1/g.V * δx⁺(i, j, k, g::grids, adv_flux_x, u, q, ΔT)
@inline div_flux_y(i, j, k, g::grids, v, q, ΔT) = 1/g.V * δy⁺(i, j, k, g::grids, adv_flux_y, v, q, ΔT)
@inline div_flux_z(i, j, k, g::grids, w, q, ΔT) = 1/g.V * δz⁺(i, j, k, g::grids, adv_flux_z, w, q, ΔT)

##### apply the tendency of by multi-dimensional advection for tracer q
@kernel function calc_q⅓!(qtemp, g::grids, u, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_x(ii, jj, kk, g::grids, u, q, ΔT)
end
@kernel function calc_q⅔!(qtemp, g::grids, v, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_y(ii, jj, kk, g::grids, v, q, ΔT)
end
@kernel function calc_q¹!(qtemp, g::grids, w, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_z(ii, jj, kk, g::grids, w, q, ΔT)
end

