#####
##### update tracer fields
#####

##### apply tendency to tracer field
@kernel function apply_tendency_kernel!(tracer, Gc, consume, g::AbstractGrid)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds tracer[ii, jj, kk] += Gc[ii, jj, kk] + consume[ii, jj, kk] /volume(ii, jj, kk, g)
end
function apply_tendency!(tracers, Gcs, consume, g::AbstractGrid, arch::Architecture)
    kernel! = apply_tendency_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(tracers[name].data, Gcs[name].data, consume[name].data, g)
    end
    return nothing
end
    
function tracer_update!(tracers, Gcs, tracer_temp, arch::Architecture, g::AbstractGrid, params, vel, consume, ΔT, iter)
    ##### compute advection tendency
    tracer_advection!(tracers, tracer_temp, Gcs, vel, g, ΔT, arch)

    ##### compute tracer diffusion,for each time step
    tracer_diffusion!(Gcs, arch, g, tracers, params["κh"], params["κh"], params["κv"], ΔT)

    ##### compute biogeochemical forcings of tracers,for each time step
    zero_fields!(tracer_temp)
    tracer_forcing!(Gcs, tracer_temp, tracers, params, ΔT)

    ##### apply boundary conditions
    apply_bcs!(Gcs, tracers, g, iter, ΔT, arch)

    ##### apply diffusion and forcing tendency
    apply_tendency!(tracers, Gcs, consume, g, arch)

    fill_halo_tracer!(tracers, g)
end