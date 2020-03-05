"""
    sum_consume(consume_p, consume_z)
sum up consumptions from zooplankton update function and phytoplankton
update function
"""
function sum_consume(consume_p, consume_z)
    consume_p.DOC = consume_z.DOC .+ consume_p.DOC
    consume_p.POC = consume_z.POC .+ consume_p.POC
    consume_p.DON = consume_z.DON .+ consume_p.DON
    consume_p.PON = consume_z.PON .+ consume_p.PON
    consume_p.DOP = consume_z.DOP .+ consume_p.DOP
    consume_p.POP = consume_z.POP .+ consume_p.POP
    consume_p.DIC = consume_z.DIC .+ consume_p.DIC
    consume_p.NH4 = consume_z.NH4 .+ consume_p.NH4
    consume_p.NO3 = consume_z.NO3 .+ consume_p.NO3
    consume_p.PO4 = consume_z.PO4 .+ consume_p.PO4
    return consume_p
end

"""
    function nut_update(model, velᵇ, consume, ΔT, Dim)
Update nutrient fields to next time step with source term and consumption by phytoplankton ('consume')
"""
function nut_update(model, velᵇ, consume, ΔT, Dim=3)
    nutrients = model.nutrients
    params = model.params
    g = model.grid
    κh = params["κh"]
    κv = params["κv"]
    # compute Forcings (remineralization)
    DIC = copy(nutrients.DIC); NH4 = copy(nutrients.NH4);
    NO3 = copy(nutrients.NO3); PO4 = copy(nutrients.PO4);
    DOC = copy(nutrients.DOC); DON = copy(nutrients.DON);
    DOP = copy(nutrients.DOP); POC = copy(nutrients.POC);
    PON = copy(nutrients.PON); POP = copy(nutrients.POP);
    DIC[DIC .< 0.0] .= 0.0; NH4[NH4 .< 0.0] .= 0.0;
    NO3[NO3 .< 0.0] .= 0.0; PO4[PO4 .< 0.0] .= 0.0;
    DOC[DOC .< 0.0] .= 0.0; DON[DON .< 0.0] .= 0.0;
    DOP[DOP .< 0.0] .= 0.0; POC[POC .< 0.0] .= 0.0;
    PON[PON .< 0.0] .= 0.0; POP[POP .< 0.0] .= 0.0;
    F = nutrients_init(g)
    # compute remineralization of organic nutrients
    F.DIC .= F.DIC .+ DOC .* params["kDOC"]
    F.NH4 .= F.NH4 .+ DON .* params["kDON"] .- NH4 .* params["Nit"]
    F.NO3 .= F.NO3 .+ NH4 .* params["Nit"]
    F.PO4 .= F.PO4 .+ DOP .* params["kDOP"]
    F.DOC .= F.DOC .- DOC .* params["kDOC"] .+ POC .* params["kPOC"]
    F.DON .= F.DON .- DON .* params["kDON"] .+ PON .* params["kPON"]
    F.DOP .= F.DOP .- DOP .* params["kDOP"] .+ POP .* params["kPOP"]
    F.POC .= F.POC .- POC .* params["kPOC"]
    F.PON .= F.PON .- PON .* params["kPON"]
    F.POP .= F.POP .- POP .* params["kPOP"]

    # Compute nutrient advection using DST3FL scheme 'MultiDim_adv(g, DIC, vel, ΔT)
    # Compute nutrient diffusion 'κ∇²(g, DIC, κh, κv, i, j, k)
    gtr = nutrients_init(g)
    if Dim == 0
        nothing
    else
        adv_DIC = MultiDim_adv(g, DIC, velᵇ, ΔT);
        adv_NH4 = MultiDim_adv(g, NH4, velᵇ, ΔT);
        adv_NO3 = MultiDim_adv(g, NO3, velᵇ, ΔT);
        adv_PO4 = MultiDim_adv(g, PO4, velᵇ, ΔT);
        adv_DOC = MultiDim_adv(g, DOC, velᵇ, ΔT);
        adv_DON = MultiDim_adv(g, DON, velᵇ, ΔT);
        adv_DOP = MultiDim_adv(g, DOP, velᵇ, ΔT);
        adv_POC = MultiDim_adv(g, POC, velᵇ, ΔT);
        adv_PON = MultiDim_adv(g, PON, velᵇ, ΔT);
        adv_POP = MultiDim_adv(g, POP, velᵇ, ΔT);
        for k in 1:g.Nz
            for j in 1:g.Ny
                for i in 1:g.Nx
                    gtr.DIC[i, j, k] = adv_DIC[i, j, k] + κ∇²(g, DIC, κh, κv, i, j, k) + F.DIC[i, j, k]
                    gtr.NH4[i, j, k] = adv_NH4[i, j, k] + κ∇²(g, NH4, κh, κv, i, j, k) + F.NH4[i, j, k]
                    gtr.NO3[i, j, k] = adv_NO3[i, j, k] + κ∇²(g, NO3, κh, κv, i, j, k) + F.NO3[i, j, k]
                    gtr.PO4[i, j, k] = adv_PO4[i, j, k] + κ∇²(g, PO4, κh, κv, i, j, k) + F.PO4[i, j, k]
                    gtr.DOC[i, j, k] = adv_DOC[i, j, k] + κ∇²(g, DOC, κh, κv, i, j, k) + F.DOC[i, j, k]
                    gtr.DON[i, j, k] = adv_DON[i, j, k] + κ∇²(g, DON, κh, κv, i, j, k) + F.DON[i, j, k]
                    gtr.DOP[i, j, k] = adv_DOP[i, j, k] + κ∇²(g, DOP, κh, κv, i, j, k) + F.DOP[i, j, k]
                    gtr.POC[i, j, k] = adv_POC[i, j, k] + κ∇²(g, POC, κh, κv, i, j, k) + F.POC[i, j, k]
                    gtr.PON[i, j, k] = adv_PON[i, j, k] + κ∇²(g, PON, κh, κv, i, j, k) + F.PON[i, j, k]
                    gtr.POP[i, j, k] = adv_POP[i, j, k] + κ∇²(g, POP, κh, κv, i, j, k) + F.POP[i, j, k]
                end
            end
        end
    end
    # store nutrients of the former time step
    nutₜ = nutrients_init(g)
    nutₜ.DIC .= nutrients.DIC .+ gtr.DIC .* ΔT .+ consume.DIC ./ g.V
    nutₜ.NH4 .= nutrients.NH4 .+ gtr.DIC .* ΔT .+ consume.NH4 ./ g.V
    nutₜ.NO3 .= nutrients.NO3 .+ gtr.DIC .* ΔT .+ consume.NO3 ./ g.V
    nutₜ.PO4 .= nutrients.PO4 .+ gtr.DIC .* ΔT .+ consume.PO4 ./ g.V
    nutₜ.DOC .= nutrients.DOC .+ gtr.DOC .* ΔT .+ consume.DOC ./ g.V
    nutₜ.DON .= nutrients.DON .+ gtr.DON .* ΔT .+ consume.DON ./ g.V
    nutₜ.DOP .= nutrients.DOP .+ gtr.DOP .* ΔT .+ consume.DOP ./ g.V
    nutₜ.POC .= nutrients.POC .+ gtr.POC .* ΔT .+ consume.POC ./ g.V
    nutₜ.PON .= nutrients.PON .+ gtr.PON .* ΔT .+ consume.PON ./ g.V
    nutₜ.POP .= nutrients.POP .+ gtr.POP .* ΔT .+ consume.POP ./ g.V
    return nutₜ, gtr
end

