##### update nutrient fields
function nut_update!(nutrients, Gcs, MD1, MD2, MD3, arch::Architecture, g::Grids, params, vel, consume, ΔT)
    # compute advection tendency
    nut_advection!(Gcs, arch, g, nutrients, MD1, MD2, MD3, vel, ΔT)

    # compute nutrient diffusion,for each time step
    nut_diffusion!(Gcs, arch, g, nutrients, params["κh"], params["κh"], params["κv"], ΔT)

    # compute biogeochemical forcings of nutrients,for each time step
    nut_forcing!(Gcs, arch, g, nutrients, params, ΔT)

    # apply diffusion and forcing tendency
    for name in nut_names
        nutrients[name].data .= nutrients[name].data .+ Gcs[name].data .+ consume[name].data ./ g.V
    end

    fill_halo!(nutrients, g)
end

function nut_update!(nutrients, Gcs, arch::Architecture, g::Grids, params, consume, ΔT)
    # compute biogeochemical forcings of nutrients,for each time step
    nut_forcing!(Gcs, arch, g, nutrients, params, ΔT)

    # apply forcing tendency
    for name in nut_names
        nutrients[name].data .= nutrients[name].data .+ Gcs[name].data .+ consume[name].data ./ g.V
    end

    fill_halo!(nutrients, g)
end
