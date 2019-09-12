###################################################################################
# Compute advective flux of tracers using 3rd order DST Scheme with flux limiting #
# Due to free surface, tracers will not be conserved at each time step.           #
# Surface flux(τⁿ[:,:,1]*Az*wFld[:,:,1]) need to be recorded to check tracer      # 
# conservation.                                                                   #
# There will still be some tiny negative values in tracer field because of multi- #
# dimensional advection.                                                          #
###################################################################################
const θmax = 1.0e20
# Increment and decrement integer a with periodic wrapping
incmod1(a, n) = ifelse(a==n, 1, a+1)
decmod1(a, n) = ifelse(a==1, n, a-1)
function decmod2(a, n)
    if a == 1
        a = max(1, n - 1)
    elseif a == 2
        a = n 
    else
        a = a - 2
    end
    return a
end

# calculate CFL number: c=uΔt/Lx
calUCFL(g::grids, ΔT, uFld, i, j, k) = abs(uFld[i, j, k] * ΔT / g.dxC[i, j])
calVCFL(g::grids, ΔT, vFld, i, j, k) = abs(vFld[i, j, k] * ΔT / g.dyC[i, j])
calWCFL(g::grids, ΔT, wFld, i, j, k) = abs(wFld[i, j, k] * ΔT / g.dzC[k])

# calculate d₀ and d₁ 
d0(CFL) = (2.0 - CFL) * (1.0 - CFL) / 6.0 
d1(CFL) = (1.0 - CFL * CFL) / 6.0 

# calculate volume transport, unit: m³/s
calUTrans(g::grids, uFld, i, j, k) = g.Ax[i, j, k] * uFld[i, j, k]
calVTrans(g::grids, vFld, i, j, k) = g.Ay[i, j, k] * vFld[i, j, k]
calWTrans(g::grids, wFld, i, j, k) = g.Az[i, j] * wFld[i, j, k]
δuTrans(g::grids, uFld, i, j, k) = calUTrans(g, uFld, incmod1(i, g.Nx), j, k) - calUTrans(g, uFld, i, j, k)
δvTrans(g::grids, vFld, i, j, k) = calVTrans(g, vFld, i, incmod1(j, g.Ny), k) - calVTrans(g, vFld, i, j, k)
function δwTrans(g::grids, wFld, i, j, k)
    if k == g.Nz
        return calWTrans(g, wFld, i, j, k)
    else
        return (calWTrans(g, wFld, i, j, k) - calWTrans(g, wFld, i, j, min(k+1, g.Nz)))
    end
end

# calculate Zonal advection flux with Sweby limiter, unit: mmolC/s
function adv_x(g::grids, q,  uFld, i, j, k, ΔT)
    uCFL = calUCFL(g, ΔT, uFld, i, j, k)
    uTrans = calUTrans(g, uFld, i, j, k)
    d₀ = d0(uCFL); d₁ = d1(uCFL);
    rj⁺= q[incmod1(i, g.Nx), j, k] - q[i, j, k]
    rj = q[i, j, k] - q[decmod1(i, g.Nx), j, k]
    rj⁻= q[decmod1(i, g.Nx), j, k] - q[decmod2(i, g.Nx), j, k] 
    if abs(rj)*θmax ≤ abs(rj⁻)
        θ⁺ = copysign(θmax, rj⁻*rj)
    else
        θ⁺ = rj⁻ / rj
    end
    if abs(rj)*θmax ≤ abs(rj⁺)
        θ⁻ = copysign(θmax, rj⁺*rj)
    else
        θ⁻ = rj⁺ / rj
    end
    Ψ⁺ = d₀ + d₁ * θ⁺
    Ψ⁺ = max(0.0, min(min(1.0, Ψ⁺), θ⁺ * (1.0 - uCFL) / (uCFL + 1.0e-20)))
    Ψ⁻ = d₀ + d₁ * θ⁻
    Ψ⁻ = max(0.0, min(min(1.0, Ψ⁻), θ⁻ * (1.0 - uCFL) / (uCFL + 1.0e-20)))
    Fᵤ  = 0.5 * (uTrans + abs(uTrans))*(q[decmod1(i, g.Nx), j, k] + Ψ⁺ * rj) + 0.5 * (uTrans - abs(uTrans)) * (q[i, j, k] - Ψ⁻ * rj)
    return Fᵤ
end

# calculate Meridional advection flux with Sweby limiter, unit: mmolC/s
function adv_y(g::grids, q,  vFld, i, j, k, ΔT)
    vCFL = calVCFL(g, ΔT, vFld, i, j, k)
    vTrans = calVTrans(g, vFld, i, j, k)
    d₀ = d0(vCFL); d₁ = d1(vCFL);
    rj⁺= q[i, incmod1(j, g.Ny), k] - q[i, j, k]
    rj = q[i, j, k] - q[i, decmod1(j, g.Ny), k]
    rj⁻= q[i, decmod1(j, g.Ny), k] - q[i, decmod2(j, g.Ny), k] 
    if abs(rj)*θmax ≤ abs(rj⁻)
        θ⁺ = copysign(θmax, rj⁻*rj)
    else
        θ⁺ = rj⁻ / rj
    end
    if abs(rj)*θmax ≤ abs(rj⁺)
        θ⁻ = copysign(θmax, rj⁺*rj)
    else
        θ⁻ = rj⁺ / rj
    end
    Ψ⁺ = d₀ + d₁ * θ⁺
    Ψ⁺ = max(0.0, min(min(1.0, Ψ⁺), θ⁺ * (1.0 - vCFL) / (vCFL + 1.0e-20)))
    Ψ⁻ = d₀ + d₁ * θ⁻
    Ψ⁻ = max(0.0, min(min(1.0, Ψ⁻), θ⁻ * (1.0 - vCFL) / (vCFL + 1.0e-20)))
    Fᵥ = 0.5 * (vTrans + abs(vTrans))*(q[i, decmod1(j, g.Ny), k] + Ψ⁺ * rj) + 0.5 * (vTrans - abs(vTrans)) * (q[i, j, k] - Ψ⁻ * rj)
    return Fᵥ
end

# calculate vertical flux with Sweby limiter, unit: mmolC/s
function adv_z(g::grids, q,  wFld, i, j, k, ΔT)
    wCFL = calWCFL(g, ΔT, wFld, i, j, k)
    wTrans = calWTrans(g, wFld, i, j, k)
    d₀ = d0(wCFL); d₁ = d1(wCFL);
    km1 = max(1, k-1); km2 = max(1, k-2); kp1 = min(g.Nz, k+1)
    rj⁺= q[i, j, k] - q[i, j, kp1]
    rj = q[i, j, km1] - q[i, j, k]
    rj⁻= q[i, j, km2] - q[i, j, km1] 
    if abs(rj)*θmax ≤ abs(rj⁻)
        θ⁺ = copysign(θmax, rj⁻*rj)
    else
        θ⁺ = rj⁻ / rj
    end
    if abs(rj)*θmax ≤ abs(rj⁺)
        θ⁻ = copysign(θmax, rj⁺*rj)
    else
        θ⁻ = rj⁺ / rj
    end
    Ψ⁺ = d₀ + d₁ * θ⁺
    Ψ⁺ = max(0.0, min(min(1.0, Ψ⁺), θ⁺ * (1.0 - wCFL) / (wCFL + 1.0e-20)))
    Ψ⁻ = d₀ + d₁ * θ⁻
    Ψ⁻ = max(0.0, min(min(1.0, Ψ⁻), θ⁻ * (1.0 - wCFL) / (wCFL + 1.0e-20)))
    Fᵣ = 0.5 * (wTrans + abs(wTrans))*(q[i, j, k] + Ψ⁻ * rj) + 0.5 * (wTrans - abs(wTrans)) * (q[i, j, km1] - Ψ⁺ * rj)
    return Fᵣ
end

# multi-dimensional advetion
function MultiDim_adv(g::grids, q, vel, ΔT)
    u = vel.u; v = vel.v; w = vel.w;
    adv = zeros(size(q)); q₁ = zeros(size(q)); q₂ = zeros(size(q)); q₃ = zeros(size(q));
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                q₁[i, j, k] =  q[i, j, k] - ΔT / g.V[i, j, k] * (adv_x(g, q, u, incmod1(i, g.Nx), j, k, ΔT) - adv_x(g, q, u, i, j, k, ΔT) - q[i, j, k] * δuTrans(g, u, i, j, k))
            end
        end
    end
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                q₂[i, j, k] = q₁[i, j, k] - ΔT / g.V[i, j, k] * (adv_y(g, q₁, v, i, incmod1(j, g.Ny), k, ΔT) - adv_y(g, q₁, v, i, j, k, ΔT) - q[i, j, k] * δvTrans(g, v, i, j, k)) 
            end
        end
    end
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                if k == g.Nz
                    q₃[i, j, k] = q₂[i, j, k] - ΔT / g.V[i, j, k] * (adv_z(g, q₂, w, i, j, k, ΔT) - q[i, j, k] * δwTrans(g, w, i, j, k))
                else
                    q₃[i, j, k] = q₂[i, j, k] - ΔT / g.V[i, j, k] * (adv_z(g, q₂, w, i, j, k, ΔT) - adv_z(g, q₂, w, i, j, min(k+1, g.Nz), ΔT) - q[i, j, k] * δwTrans(g, w, i, j, k))
                end
            end
        end
    end
    adv .= (q₃ .- q) ./ ΔT
    return adv
end
