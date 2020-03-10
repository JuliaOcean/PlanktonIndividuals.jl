"""
    sum_nut_tendency(consume_p, consume_z)
sum up 2 nut tendencies
"""
function sum_nut_tendency(a::nutrient_fields, b::nutrient_fields)
    a.DOC = b.DOC .+ a.DOC
    a.POC = b.POC .+ a.POC
    a.DON = b.DON .+ a.DON
    a.PON = b.PON .+ a.PON
    a.DOP = b.DOP .+ a.DOP
    a.POP = b.POP .+ a.POP
    a.DIC = b.DIC .+ a.DIC
    a.NH4 = b.NH4 .+ a.NH4
    a.NO3 = b.NO3 .+ a.NO3
    a.PO4 = b.PO4 .+ a.PO4
    return a
end

"""
    nut_forcing(nut, g, params)
compute remineralization and nitrification of nutrients
"""
function nut_forcing(g, nutrients, params,ΔT)
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
    F.DIC .= F.DIC .+ DOC .* params["kDOC"] .* ΔT
    F.DOC .= F.DOC .- DOC .* params["kDOC"] .* ΔT
    F.DOC .= F.DOC .+ POC .* params["kPOC"] .* ΔT
    F.POC .= F.POC .- POC .* params["kPOC"] .* ΔT

    F.NH4 .= F.NH4 .+ DON .* params["kDON"] .* ΔT
    F.NH4 .= F.NH4 .- NH4 .* params["Nit"]  .* ΔT
    F.NO3 .= F.NO3 .+ NH4 .* params["Nit"]  .* ΔT
    F.DON .= F.DON .- DON .* params["kDON"] .* ΔT
    F.DON .= F.DON .+ PON .* params["kPON"] .* ΔT
    F.PON .= F.PON .- PON .* params["kPON"] .* ΔT

    F.PO4 .= F.PO4 .+ DOP .* params["kDOP"] .* ΔT
    F.DOP .= F.DOP .- DOP .* params["kDOP"] .* ΔT
    F.DOP .= F.DOP .+ POP .* params["kPOP"] .* ΔT
    F.POP .= F.POP .- POP .* params["kPOP"] .* ΔT
    return F
end

"""
    nut_advection(g, nutrients, velᵇ, ΔT)
compute nutrient advection using DST3FL shceme
"""
function nut_advection(g, nutrients, velᵇ, ΔT)
    gtr = nutrients_init(g)
    gtr.DIC = MultiDim_adv(g, nutrients.DIC, velᵇ, ΔT) .* ΔT;
    gtr.NH4 = MultiDim_adv(g, nutrients.NH4, velᵇ, ΔT) .* ΔT;
    gtr.NO3 = MultiDim_adv(g, nutrients.NO3, velᵇ, ΔT) .* ΔT;
    gtr.PO4 = MultiDim_adv(g, nutrients.PO4, velᵇ, ΔT) .* ΔT;
    gtr.DOC = MultiDim_adv(g, nutrients.DOC, velᵇ, ΔT) .* ΔT;
    gtr.DON = MultiDim_adv(g, nutrients.DON, velᵇ, ΔT) .* ΔT;
    gtr.DOP = MultiDim_adv(g, nutrients.DOP, velᵇ, ΔT) .* ΔT;
    gtr.POC = MultiDim_adv(g, nutrients.POC, velᵇ, ΔT) .* ΔT;
    gtr.PON = MultiDim_adv(g, nutrients.PON, velᵇ, ΔT) .* ΔT;
    gtr.POP = MultiDim_adv(g, nutrients.POP, velᵇ, ΔT) .* ΔT;
    return gtr
end

"""
    nut_diffusion(g, nutrients, params)
compute diffusion for each nutrient tracer
"""
function nut_diffusion(g, nutrients, params, ΔT)
    diffu = nutrients_init(g)
    κh = params["κh"]
    κv = params["κv"]
    for k in 1:g.Nz
        for j in 1:g.Ny
            for i in 1:g.Nx
                diffu.DIC[i, j, k] =κ∇²(g, nutrients.DIC, κh, κv, i, j, k) * ΔT
                diffu.NH4[i, j, k] =κ∇²(g, nutrients.NH4, κh, κv, i, j, k) * ΔT
                diffu.NO3[i, j, k] =κ∇²(g, nutrients.NO3, κh, κv, i, j, k) * ΔT
                diffu.PO4[i, j, k] =κ∇²(g, nutrients.PO4, κh, κv, i, j, k) * ΔT
                diffu.DOC[i, j, k] =κ∇²(g, nutrients.DOC, κh, κv, i, j, k) * ΔT
                diffu.DON[i, j, k] =κ∇²(g, nutrients.DON, κh, κv, i, j, k) * ΔT
                diffu.DOP[i, j, k] =κ∇²(g, nutrients.DOP, κh, κv, i, j, k) * ΔT
                diffu.POC[i, j, k] =κ∇²(g, nutrients.POC, κh, κv, i, j, k) * ΔT
                diffu.PON[i, j, k] =κ∇²(g, nutrients.PON, κh, κv, i, j, k) * ΔT
                diffu.POP[i, j, k] =κ∇²(g, nutrients.POP, κh, κv, i, j, k) * ΔT
            end
        end
    end
    return diffu
end

"""
    function nut_update(model, velᵇ, consume, ΔT, Dim)
Update nutrient fields to next time step with source term and consumption by phytoplankton ('consume')
"""
function nut_update(model, velᵇ, consume, ΔT)
    nutrients = model.nutrients
    params = model.params
    g = model.grid
    # compute biogeochemical forcings of nutrients,for each time step
    F = nut_forcing(g, nutrients, params, ΔT)
    # Compute nutrient advection using DST3FL scheme,for each time step
    gtr = nut_advection(g, nutrients, velᵇ, ΔT)
    # Compute nutrient diffusion,for each time step
    diffu = nut_diffusion(g, nutrients, params, ΔT)
    # sum all tendencies
    gtr = sum_nut_tendency(gtr,diffu)
    tendencies = sum_nut_tendency(gtr,F)
    # tendencies = gtr
    # store nutrients of the former time step
    nutₜ = nutrients_init(g)
    nutₜ.DIC .= nutrients.DIC .+ tendencies.DIC .+ consume.DIC ./ g.V
    nutₜ.NH4 .= nutrients.NH4 .+ tendencies.NH4 .+ consume.NH4 ./ g.V
    nutₜ.NO3 .= nutrients.NO3 .+ tendencies.NO3 .+ consume.NO3 ./ g.V
    nutₜ.PO4 .= nutrients.PO4 .+ tendencies.PO4 .+ consume.PO4 ./ g.V
    nutₜ.DOC .= nutrients.DOC .+ tendencies.DOC .+ consume.DOC ./ g.V
    nutₜ.DON .= nutrients.DON .+ tendencies.DON .+ consume.DON ./ g.V
    nutₜ.DOP .= nutrients.DOP .+ tendencies.DOP .+ consume.DOP ./ g.V
    nutₜ.POC .= nutrients.POC .+ tendencies.POC .+ consume.POC ./ g.V
    nutₜ.PON .= nutrients.PON .+ tendencies.PON .+ consume.PON ./ g.V
    nutₜ.POP .= nutrients.POP .+ tendencies.POP .+ consume.POP ./ g.V
    return nutₜ, gtr
end

function nut_update(model, consume, ΔT)
    nutrients = model.nutrients
    params = model.params
    g = model.grid
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

    gtr = nutrients_init(g)
    # store nutrients of the former time step
    nutₜ = nutrients_init(g)
    nutₜ.DIC .= nutrients.DIC .+ gtr.DIC .* ΔT .+ consume.DIC ./ g.V
    nutₜ.NH4 .= nutrients.NH4 .+ gtr.NH4 .* ΔT .+ consume.NH4 ./ g.V
    nutₜ.NO3 .= nutrients.NO3 .+ gtr.NO3 .* ΔT .+ consume.NO3 ./ g.V
    nutₜ.PO4 .= nutrients.PO4 .+ gtr.PO4 .* ΔT .+ consume.PO4 ./ g.V
    nutₜ.DOC .= nutrients.DOC .+ gtr.DOC .* ΔT .+ consume.DOC ./ g.V
    nutₜ.DON .= nutrients.DON .+ gtr.DON .* ΔT .+ consume.DON ./ g.V
    nutₜ.DOP .= nutrients.DOP .+ gtr.DOP .* ΔT .+ consume.DOP ./ g.V
    nutₜ.POC .= nutrients.POC .+ gtr.POC .* ΔT .+ consume.POC ./ g.V
    nutₜ.PON .= nutrients.PON .+ gtr.PON .* ΔT .+ consume.PON ./ g.V
    nutₜ.POP .= nutrients.POP .+ gtr.POP .* ΔT .+ consume.POP ./ g.V
    return nutₜ, gtr
end
