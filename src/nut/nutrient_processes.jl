"""
    function nut_update(model, velᵇ, consume, ΔT, Dim)
Update nutrient fields to next time step with source term and consumption by phytoplankton ('consume')
"""
function nut_update_interior(model, velᵇ, consume, ΔT)
    nutrients = model.nutrients
    params = model.params
    g = model.grid
    arch = model.arch

    # compute biogeochemical forcings of nutrients,for each time step
    F = nutrients_init(arch, g)
    nut_forcing!(F, arch, g, nutrients, params, ΔT)

    # compute nutrient diffusion,for each time step
    diffu = nutrients_init(arch, g)
    nut_diffusion!(diffu, arch, g, nutrients, params["κh"], params["κh"], params["κv"], ΔT)

    # compute advection tendency
    nutₜ = nutrients_init(arch, g)
    nut_advection!(nutₜ, arch, g, nutrients, velᵇ, ΔT)

    # compute advection tendency for record
    gtr = nutrients_init(arch, g)
    gtr = sub_nut_tendency!(gtr, nutₜ,nutrients)

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

    fill_halo!(nutₜ, g)

    return nutₜ, add_nut_tendency!(gtr, add_nut_tendency!(diffu, F))
end

function nut_update_interior(model, consume, ΔT)
    nutrients = model.nutrients
    params = model.params
    g = model.grid
    arch = model.arch

    # compute biogeochemical forcings of nutrients,for each time step
    F = nutrients_init(arch, g)
    nut_forcing!(F, arch, g, nutrients, params, ΔT)

    # apply forcing tendency
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

    fill_halo!(nutₜ, g)

    return nutₜ, F
end
