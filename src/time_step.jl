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
    agent_num = size(phyts_b,1)
    write_pop_dynamics(model.t, agent_num, counts_p, resultspath)
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
    model.individuals.phytos = phyts_b
    agent_num = size(phyts_b,1)
    write_pop_dynamics(model.t, agent_num, counts_p, resultspath)
    nutₜ,gtr = nut_update(model, consume_p, ΔT)
    write_nut_cons(model.grid, gtr, nutₜ,model.t,resultspath)
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
    for j in 1:size(planks,1)
        plank = planks[j,:]
        agent_advection(plank,vel_itp,grid,ΔT)
        if grid.Nx > 1
            agent_diffusionX(plank,grid,params["κhP"])
        end
        if grid.Ny > 1
            agent_diffusionY(plank,grid,params["κhP"])
        end
        if grid.Nz > 1
            agent_diffusionZ(plank,grid,params["κvP"])
        end
    end
    if model.individuals.zoos ≠ nothing
        zoos = model.individuals.zoos
        for j in 1:size(zoos,1)
            plank = zoos[j,:]
            agent_advection(plank,vel_itp,grid,ΔT)
            if grid.Nx > 1
                agent_diffusionX(plank,grid,params["κhP"])
            end
            if grid.Ny > 1
                agent_diffusionY(plank,grid,params["κhP"])
            end
            if grid.Nz > 1
                agent_diffusionZ(plank,grid,params["κvP"])
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
    for j in 1:size(planks,1)
        plank = planks[j,:]
        agent_advectionRK4(plank,vel_itps,grid,ΔT)
        if grid.Nx > 1
            agent_diffusionX(plank,grid,params["κhP"])
        end
        if grid.Ny > 1
            agent_diffusionY(plank,grid,params["κhP"])
        end
        if grid.Nz > 1
            agent_diffusionZ(plank,grid,params["κvP"])
        end
    end
    if model.individuals.zoos ≠ nothing
        zoos = model.individuals.zoos
        for j in 1:size(zoos,1)
            plank = zoos[j,:]
            agent_advectionRK4(plank,vel_itps,grid,ΔT)
            if grid.Nx > 1
                agent_diffusionX(plank,grid,params["κhP"])
            end
            if grid.Ny > 1
                agent_diffusionY(plank,grid,params["κhP"])
            end
            if grid.Nz > 1
                agent_diffusionZ(plank,grid,params["κvP"])
            end
        end
    end
end
