#####
##### update nutrient fields
#####

##### apply tendency to nutrient field
@kernel function apply_tendency_kernel!(nut, Gc, consume, ΔT, g::AbstractGrid)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds nut[ii, jj, kk] += Gc[ii, jj, kk] + consume[ii, jj, kk] /volume(ii, jj, kk, g)
end
function apply_tendency!(nuts, Gcs, consume, ΔT, g::AbstractGrid, arch::Architecture)
    kernel! = apply_tendency_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nuts[name].data, Gcs[name].data, consume[name].data, ΔT, g, dependencies=barrier)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end
    
function nut_update!(nutrients, Gcs, nut_temp, arch::Architecture, g::AbstractGrid, params, vel, consume, ΔT, iter)
    ##### compute advection tendency
    nut_advection!(nutrients, nut_temp, Gcs, vel, g, ΔT, arch)

    ##### compute nutrient diffusion,for each time step
    nut_diffusion!(Gcs, arch, g, nutrients, params["κh"], ΔT)

    ##### compute biogeochemical forcings of nutrients,for each time step
    zero_fields!(nut_temp)
    nut_forcing!(Gcs, nut_temp, nutrients, params, ΔT)

    ##### apply boundary conditions
    apply_bcs!(Gcs, nutrients, g, iter, arch)

    ##### apply diffusion and forcing tendency
    apply_tendency!(nutrients, Gcs, consume, ΔT, g, arch)

    fill_halo_nut!(nutrients, g)
end

function nut_update!(nutrients, Gcs, nut_temp, arch::Architecture, g::AbstractGrid, params, consume, ΔT, iter)
    # compute biogeochemical forcings of nutrients,for each time step
    nut_forcing!(Gcs, nut_temp, nutrients, params, ΔT)

    ##### apply boundary conditions
    apply_bcs!(Gcs, nutrients, g, iter, arch)

    # apply forcing tendency
    apply_tendency!(nutrients, Gcs, consume, ΔT, g, arch)

    fill_halo_nut!(nutrients, g)
end
