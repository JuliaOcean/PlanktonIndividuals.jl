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
    add_nut_tendency!(diffu, F)

    # compute advection tendency
    gtr = nutrients_init(arch, g)
    nut_advection!(gtr, arch, g, nutrients, velᵇ, ΔT)
    add_nut_tendency!(gtr, diffu)

    # apply diffusion and forcing tendency
    nutₜ = nutrients_init(arch, g)
    for name in nut_names
        nutₜ[name].data .= nutrients[name].data .+ gtr[name].data .+ consume[name].data ./ g.V
    end

    fill_halo!(nutₜ, g)

    return nutₜ, gtr
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

    for name in nut_names
        nutₜ[name].data .= nutrients[name].data .+ F[name].data .+ consume[name].data ./ g.V
    end

    fill_halo!(nutₜ, g)

    return nutₜ, F
end
