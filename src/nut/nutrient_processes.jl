##### update nutrient fields
function nut_update!(nutrients, gtr, arch::Architecture, g::Grids, params, vel, consume, ΔT)
    # compute biogeochemical forcings of nutrients,for each time step
    F = nutrients_init(arch, g)
    nut_forcing!(F, arch, g, nutrients, params, ΔT)

    # compute nutrient diffusion,for each time step
    diffu = nutrients_init(arch, g)
    nut_diffusion!(diffu, arch, g, nutrients, params["κh"], params["κh"], params["κv"], ΔT)
    add_nut_tendency!(diffu, F)

    # compute advection tendency
    nut_advection!(gtr, arch, g, nutrients, vel, ΔT)
    add_nut_tendency!(gtr, diffu)

    # apply diffusion and forcing tendency
    for name in nut_names
        nutrients[name].data .= nutrients[name].data .+ gtr[name].data .+ consume[name].data ./ g.V
    end

    fill_halo!(nutrients, g)
end

function nut_update!(nutrients, gtr, arch::Architecture, g::Grids, params, consume, ΔT)
    # compute biogeochemical forcings of nutrients,for each time step
    nut_forcing!(gtr, arch, g, nutrients, params, ΔT)

    # apply forcing tendency
    for name in nut_names
        nutrients[name].data .= nutrients[name].data .+ gtr[name].data .+ consume[name].data ./ g.V
    end

    fill_halo!(nutrients, g)
end
