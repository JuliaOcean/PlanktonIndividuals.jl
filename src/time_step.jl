"""
    PI_TimeStep!(model, RunParam, velᵇ)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PI_TimeStep!(model::Model_struct, ΔT, velᵇ::velocity, resultspath)
    model.t += 1
    phyts_b,counts_p,consume_p=phyt_update(model, ΔT)
    model.individuals.phytos = phyts_b
    if model.individuals.zoos ≠ nothing
        zoos_b,counts_z,consume_z =zoo_update(model, ΔT)
        consume_p = sum_nut_tendency(consume_p, consume_z)
        model.individuals.zoos = zoos_b
        counts_p.graze = counts_z.graze
    end
    write_pop_dynamics(model.t, counts_p, resultspath)
    nutₜ,gtr = nut_update(model, velᵇ, consume_p, ΔT)
    write_nut_cons(model.grid, gtr, nutₜ,model.t,resultspath)
    model.nutrients = nutₜ
end
function PI_TimeStep!(model::Model_struct, ΔT, velᵇ::velocity)
    model.t += 1
    phyts_b,counts_p,consume_p=phyt_update(model, ΔT)
    model.individuals.phytos = phyts_b
    if model.individuals.zoos ≠ nothing
        zoos_b,counts_z,consume_z =zoo_update(model, ΔT)
        consume_p = sum_nut_tendency(consume_p, consume_z)
        model.individuals.zoos = zoos_b
        counts_p.graze = counts_z.graze
    end
    nutₜ,gtr = nut_update(model, velᵇ, consume_p, ΔT)
    model.nutrients = nutₜ
end
function PI_TimeStep!(model::Model_struct, ΔT, resultspath)
    model.t += 1
    phyts_b,counts_p,consume_p=phyt_update(model, ΔT)
    write_pop_dynamics(model.t, counts_p, resultspath)
    nutₜ,gtr = nut_update(model, consume_p, ΔT)
    write_nut_cons(model.grid, gtr, nutₜ,model.t,resultspath)
    model.individuals.phytos = phyts_b
    model.nutrients = nutₜ
end

"""
    PI_advect!(individuals, ΔT, vel)
Individual advection using a simple scheme
Particle sinking included
"""
function PI_advect!(model, ΔT, vel_itp)
    grid = model.grid
    params = model.params
    planks = model.individuals.phytos
    for j in 1:size(planks,2)
        planks[1:3,j] = agent_advection(planks[1:3,j],vel_itp,grid,ΔT)
        if grid.Nx > 1
            planks[1,j] = agent_diffusionX(planks[1,j],grid,params["κhP"])
        end
        if grid.Ny > 1
            planks[2,j] = agent_diffusionY(planks[2,j],grid,params["κhP"])
        end
        if grid.Nz > 1
            planks[3,j] = agent_diffusionZ(planks[3,j],grid,params["κvP"])
        end
    end
    if model.individuals.zoos ≠ nothing
        zoos = model.individuals.zoos
        for j in 1:size(zoos,2)
            zoos[1:3,j] = agent_advection(zoos[1:3,j],vel_itp,grid,ΔT)
            if grid.Nx > 1
                zoos[1,j] = agent_diffusionX(zoos[1,j],grid,params["κhP"])
            end
            if grid.Ny > 1
                zoos[2,j] = agent_diffusionY(zoos[2,j],grid,params["κhP"])
            end
            if grid.Nz > 1
                zoos[3,j] = agent_diffusionZ(zoos[3,j],grid,params["κvP"])
            end
        end
    end
end

"""
    PI_advectRK4!(individuals, ΔT, vel_field, grid)
Individual advection using a RK4 method
Used for 3D double grids
Particle sinking not included
"""
function PI_advectRK4!(model, ΔT, vel_itps)
    grid = model.grid
    params = model.params
    planks = model.individuals.phytos
    for j in 1:size(planks,2)
        planks[1:3,j] = agent_advectionRK4(planks[1:3,j],vel_itps,grid,ΔT)
        if grid.Nx > 1
            planks[1,j] = agent_diffusionX(planks[1,j],grid,params["κhP"])
        end
        if grid.Ny > 1
            planks[2,j] = agent_diffusionY(planks[2,j],grid,params["κhP"])
        end
        if grid.Nz > 1
            planks[3,j] = agent_diffusionZ(planks[3,j],grid,params["κvP"])
        end
    end
    if model.individuals.zoos ≠ nothing
        zoos = model.individuals.zoos
        for j in 1:size(zoos,2)
            zoos[1:3,j] = agent_advectionRK4(zoos[1:3,j],vel_itps,grid,ΔT)
            if grid.Nx > 1
                zoos[1,j] = agent_diffusionX(zoos[1,j],grid,params["κhP"])
            end
            if grid.Ny > 1
                zoos[2,j] = agent_diffusionY(zoos[2,j],grid,params["κhP"])
            end
            if grid.Nz > 1
                zoos[3,j] = agent_diffusionZ(zoos[3,j],grid,params["κvP"])
            end
        end
    end
end
