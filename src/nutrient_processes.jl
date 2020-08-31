"""
    sub_nut_tendency!(consume_p, consume_z)
subtract one tendency from total tendencies
"""
function sub_nut_tendency!(a::nutrient_fields, b::nutrient_fields, c::nutrient_fields)
    a.DOC = b.DOC .- c.DOC
    a.POC = b.POC .- c.POC
    a.DON = b.DON .- c.DON
    a.PON = b.PON .- c.PON
    a.DOP = b.DOP .- c.DOP
    a.POP = b.POP .- c.POP
    a.DIC = b.DIC .- c.DIC
    a.NH4 = b.NH4 .- c.NH4
    a.NO3 = b.NO3 .- c.NO3
    a.PO4 = b.PO4 .- c.PO4
end

"""
    nut_forcing(nut, g, params)
compute remineralization and nitrification of nutrients
"""
function nut_forcing(arch::Architecture, g, nutrients, params,ΔT)
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
    F = nutrients_init(arch, g)
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
    function nut_update(model, velᵇ, consume, ΔT, Dim)
Update nutrient fields to next time step with source term and consumption by phytoplankton ('consume')
"""
function nut_update(arch::Architecture, model, velᵇ, consume, ΔT)
    nutrients = model.nutrients
    params = model.params
    g = model.grid

    # compute biogeochemical forcings of nutrients,for each time step
    F = nut_forcing(arch, g, nutrients, params, ΔT)

    # compute nutrient diffusion,for each time step
    diffu = nutrients_init(arch, g)
    nut_diffusion!(diffu, arch, g, nutrients, params["κh"], params["κh"], params["κv"], ΔT)

    # compute advection tendency
    nutₜ = nutrients_init(arch, g)
    nut_advection!(nutₜ, arch, g, nutrients, velᵇ, ΔT)

    # apply diffusion and forcing tendency
    nutₜ.DIC .= nutₜ.DIC .+ diffu.DIC .+ F.DIC .+ consume.DIC ./ g.V
    nutₜ.NH4 .= nutₜ.NH4 .+ diffu.NH4 .+ F.NH4 .+ consume.NH4 ./ g.V
    nutₜ.NO3 .= nutₜ.NO3 .+ diffu.NO3 .+ F.NO3 .+ consume.NO3 ./ g.V
    nutₜ.PO4 .= nutₜ.PO4 .+ diffu.PO4 .+ F.PO4 .+ consume.PO4 ./ g.V
    nutₜ.DOC .= nutₜ.DOC .+ diffu.DOC .+ F.DOC .+ consume.DOC ./ g.V
    nutₜ.DON .= nutₜ.DON .+ diffu.DON .+ F.DON .+ consume.DON ./ g.V
    nutₜ.DOP .= nutₜ.DOP .+ diffu.DOP .+ F.DOP .+ consume.DOP ./ g.V
    nutₜ.POC .= nutₜ.POC .+ diffu.POC .+ F.POC .+ consume.POC ./ g.V
    nutₜ.PON .= nutₜ.PON .+ diffu.PON .+ F.PON .+ consume.PON ./ g.V
    nutₜ.POP .= nutₜ.POP .+ diffu.POP .+ F.POP .+ consume.POP ./ g.V

    # compute total tendencies
    gtr = nutrients_init(arch, g)
    gtr = sub_nut_tendency!(gtr, nutₜ,nutrients)

    return nutₜ, gtr
end

function nut_update(arch::Architecture, model, consume, ΔT)
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
    F = nutrients_init(arch, g)
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
    # store nutrients of the former time step
    nutₜ = nutrients_init(arch, g)
    nutₜ.DIC .= nutrients.DIC .+ F.DIC .+ consume.DIC ./ g.V
    nutₜ.NH4 .= nutrients.NH4 .+ F.NH4 .+ consume.NH4 ./ g.V
    nutₜ.NO3 .= nutrients.NO3 .+ F.NO3 .+ consume.NO3 ./ g.V
    nutₜ.PO4 .= nutrients.PO4 .+ F.PO4 .+ consume.PO4 ./ g.V
    nutₜ.DOC .= nutrients.DOC .+ F.DOC .+ consume.DOC ./ g.V
    nutₜ.DON .= nutrients.DON .+ F.DON .+ consume.DON ./ g.V
    nutₜ.DOP .= nutrients.DOP .+ F.DOP .+ consume.DOP ./ g.V
    nutₜ.POC .= nutrients.POC .+ F.POC .+ consume.POC ./ g.V
    nutₜ.PON .= nutrients.PON .+ F.PON .+ consume.PON ./ g.V
    nutₜ.POP .= nutrients.POP .+ F.POP .+ consume.POP ./ g.V
    return nutₜ, F
end
