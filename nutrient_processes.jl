function compute_nut_biochem(nutrients, F, rem)
    DIC, DIN, DOC, DON, POC, PON = nutrients.DIC, nutrients.DIN, nutrients.DOC, nutrients.DON, nutrients.POC, nutrients.PON
    # compute remineralization of organic nutrients
    F.DIC .= F.DIC .+ DOC .* rem.DOC
    F.DIN .= F.DIN .+ DON .* rem.DON
    F.DOC .= F.DOC .- DOC .* rem.DOC .+ POC .* rem.POC
    F.DON .= F.DON .- DON .* rem.DON .+ PON .* rem.PON
    F.POC .= F.POC .- POC .* rem.POC
    F.PON .= F.PON .- PON .* rem.PON
end



function compute_source_term(gtr, nutrients, velᵇ, g, F)
    u, v, w = velᵇ.u, velᵇ.v, velᵇ.w
    DIC, DIN, DOC, DON, POC, PON = nutrients.DIC, nutrients.DIN, nutrients.DOC, nutrients.DON, nutrients.POC, nutrients.PON
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                gtr.DIC[j, i, k] = -div_flux(g, u, v, w, DIC, i, j, k) + κ∇²(g, DIC, κh, κv, i, j, k) + F.DIC[j, i, k]
                gtr.DIN[j, i, k] = -div_flux(g, u, v, w, DIN, i, j, k) + κ∇²(g, DIN, κh, κv, i, j, k) + F.DIN[j, i, k]
                gtr.DOC[j, i, k] = -div_flux(g, u, v, w, DOC, i, j, k) + κ∇²(g, DOC, κh, κv, i, j, k) + F.DOC[j, i, k]
                gtr.DON[j, i, k] = -div_flux(g, u, v, w, DON, i, j, k) + κ∇²(g, DON, κh, κv, i, j, k) + F.DON[j, i, k]
                gtr.POC[j, i, k] = -div_flux(g, u, v, w, POC, i, j, k) + κ∇²(g, POC, κh, κv, i, j, k) + F.POC[j, i, k]
                gtr.PON[j, i, k] = -div_flux(g, u, v, w, PON, i, j, k) + κ∇²(g, PON, κh, κv, i, j, k) + F.PON[j, i, k]
            end
        end
    end
end



function nut_update(nutrients, consume, g, gtr, ΔT)
    # store nutrients of the former time step
    nutp = nutrient_fields(zeros(g.Ny, g.Nx, g.Nz), zeros(g.Ny, g.Nx, g.Nz), zeros(g.Ny, g.Nx, g.Nz), zeros(g.Ny, g.Nx, g.Nz), zeros(g.Ny, g.Nx, g.Nz), zeros(g.Ny, g.Nx, g.Nz))
    nutp.DIC = nutrients.DIC; nutp.DIN = nutrients.DIN; nutp.DOC = nutrients.DOC;
    nutp.DON = nutrients.DON; nutp.POC = nutrients.POC; nutp.PON = nutrients.PON;
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                nutrients.DIC[j, i, k] = nutrients.DIC[j, i, k] + gtr.DIC[j, i, k] * ΔT + consume.DIC[j, i, k]
                nutrients.DIN[j, i, k] = nutrients.DIN[j, i, k] + gtr.DIN[j, i, k] * ΔT + consume.DIN[j, i, k]
                nutrients.DOC[j, i, k] = nutrients.DOC[j, i, k] + gtr.DOC[j, i, k] * ΔT + consume.DOC[j, i, k]
                nutrients.DON[j, i, k] = nutrients.DON[j, i, k] + gtr.DON[j, i, k] * ΔT + consume.DON[j, i, k]
                nutrients.POC[j, i, k] = nutrients.POC[j, i, k] + gtr.POC[j, i, k] * ΔT + consume.POC[j, i, k]
                nutrients.PON[j, i, k] = nutrients.PON[j, i, k] + gtr.PON[j, i, k] * ΔT + consume.PON[j, i, k]
            end
        end
    end
    return nutp
end

