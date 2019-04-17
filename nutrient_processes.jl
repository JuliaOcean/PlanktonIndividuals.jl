function compute_nut_biochem(nutrients, rem)
    DIC, DIN, DOC, DON, POC, PON = nutrients.DIC, nutrients.DIN, nutrients.DOC, nutrients.DON, nutrients.POC, nutrients.PON
    F = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    # compute remineralization of organic nutrients
    F.DIC .= F.DIC .+ DOC .* rem.DOC
    F.DIN .= F.DIN .+ DON .* rem.DON
    F.DOC .= F.DOC .- DOC .* rem.DOC .+ POC .* rem.POC
    F.DON .= F.DON .- DON .* rem.DON .+ PON .* rem.PON
    F.POC .= F.POC .- POC .* rem.POC
    F.PON .= F.PON .- PON .* rem.PON
    return F
end



function compute_source_term(nutrients, velᵇ, g, F)
    u, v, w = velᵇ.u, velᵇ.v, velᵇ.w
    DIC, DIN, DOC, DON, POC, PON = nutrients.DIC, nutrients.DIN, nutrients.DOC, nutrients.DON, nutrients.POC, nutrients.PON
    gtr = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                gtr.DIC[i, j, k] = -div_flux(g, u, v, w, DIC, i, j, k) + κ∇²(g, DIC, κh, κv, i, j, k) + F.DIC[i, j, k]
                gtr.DIN[i, j, k] = -div_flux(g, u, v, w, DIN, i, j, k) + κ∇²(g, DIN, κh, κv, i, j, k) + F.DIN[i, j, k]
                gtr.DOC[i, j, k] = -div_flux(g, u, v, w, DOC, i, j, k) + κ∇²(g, DOC, κh, κv, i, j, k) + F.DOC[i, j, k]
                gtr.DON[i, j, k] = -div_flux(g, u, v, w, DON, i, j, k) + κ∇²(g, DON, κh, κv, i, j, k) + F.DON[i, j, k]
                gtr.POC[i, j, k] = -div_flux(g, u, v, w, POC, i, j, k) + κ∇²(g, POC, κh, κv, i, j, k) + F.POC[i, j, k]
                gtr.PON[i, j, k] = -div_flux(g, u, v, w, PON, i, j, k) + κ∇²(g, PON, κh, κv, i, j, k) + F.PON[i, j, k]
            end
        end
    end
    return gtr
end



function nut_update(nutrients, consume, g, gtr, ΔT)
    # store nutrients of the former time step
    nutₜ = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    ε = 0.0
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                nutₜ.DIC[i, j, k] = max(ε, nutrients.DIC[i, j, k] + gtr.DIC[i, j, k] * ΔT + consume.DIC[i, j, k])
                nutₜ.DIN[i, j, k] = max(ε, nutrients.DIN[i, j, k] + gtr.DIN[i, j, k] * ΔT + consume.DIN[i, j, k])
                nutₜ.DOC[i, j, k] = max(ε, nutrients.DOC[i, j, k] + gtr.DOC[i, j, k] * ΔT + consume.DOC[i, j, k])
                nutₜ.DON[i, j, k] = max(ε, nutrients.DON[i, j, k] + gtr.DON[i, j, k] * ΔT + consume.DON[i, j, k])
                nutₜ.POC[i, j, k] = max(ε, nutrients.POC[i, j, k] + gtr.POC[i, j, k] * ΔT + consume.POC[i, j, k])
                nutₜ.PON[i, j, k] = max(ε, nutrients.PON[i, j, k] + gtr.PON[i, j, k] * ΔT + consume.PON[i, j, k])
            end
        end
    end
    return nutₜ
end

