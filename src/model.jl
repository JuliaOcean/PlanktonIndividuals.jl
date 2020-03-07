"""
    PA_model(grid, RunParam)
Generate the model structure for time step
Default distribution of individuals is Normal distribution with 1.0 as mean and 0.25 as SD
Default PAR and temp are from ../samples
"""
function PA_Model(grid, RunParam;
                         t = 0,
               individuals = setup_agents(RunParam,grid),
                 nutrients,
                       PAR = read_IR_input(RunParam.nTime, RunParam.ΔT, grid),
                      temp = read_temp_input(RunParam.nTime, RunParam.ΔT, grid),
                    params = param_default,
                   )
    model = Model_struct(t,individuals, nutrients, grid, PAR, temp, params)
    if RunParam.Zoo == true
        model.params["Grz_P"] = 0
    else
        nothing
    end
    return model
end
"""
    PA_ModelRun(model::Model_struct, Rumparam, RunOption)
The function to run the model for number of time steps in RunParam
'model' is generated from 'PA_model'
The function use 'PA_advect!', a simple advection for individuals
"""
function PA_ModelRun(model::Model_struct, RunParam::RunParams, RunOption::RunOptions)
    if RunOption.NutOutputChoice == false
        DIC = zeros(g.Nx, g.Ny, g.Nz, nTime)
        NH4 = zeros(g.Nx, g.Ny, g.Nz, nTime)
        NO3 = zeros(g.Nx, g.Ny, g.Nz, nTime)
        PO4 = zeros(g.Nx, g.Ny, g.Nz, nTime)
        DOC = zeros(g.Nx, g.Ny, g.Nz, nTime)
        DON = zeros(g.Nx, g.Ny, g.Nz, nTime)
        DOP = zeros(g.Nx, g.Ny, g.Nz, nTime)
        POC = zeros(g.Nx, g.Ny, g.Nz, nTime)
        PON = zeros(g.Nx, g.Ny, g.Nz, nTime)
        POP = zeros(g.Nx, g.Ny, g.Nz, nTime)
    end
    if RunOption.VelChoice
        store_vel=load(dirname(pathof(PhytoAgentModel))*"/../samples/uvw.jld", "uvw")
    end

    for it in 1:RunParam.nTime
        t = model.t
        phyts_a = copy(model.individuals[end]) # read data from last time step
        phyts_b,counts,consume=phyt_update(t, RunParam.ΔT, phyts_a, model)
        if RunOption.VelChoice == false
            velᵇ = read_offline_vels(RunOption.VelOfflineOpt,trunc(Int,t*RunParam.ΔT/3600))
        else
            velᵇ=store_vel[t]
        end
        if (model.grid.Nx == 1) & (model.grid.Ny == 1) & (model.grid.Nz == 1)
            nothing #for 0D only
        else
            vel_itp = generate_vel_itp(model.grid, velᵇ)
            PA_advect!(model, RunParam.ΔT, vel_itp)
        end
        push!(model.individuals,phyts_b)
        write_output(t,phyts_b,counts,model.output)
        agent_num = size(phyts_b,1)
        nutₜ,gtr = nut_update(model, velᵇ, consume, RunParam.ΔT)
        if RunOption.OutputChoice
            write_nut_cons(model.grid, gtr, nutₜ, velᵇ, agent_num, t, death_ct, graz_ct, dvid_ct)
            if RunOption.NutOutputChoice
                write_nut_nc_each_step(model.grid, nutₜ, t)
            else
                DIC[:,:,:,t] = nutₜ.DIC
                NH4[:,:,:,t] = nutₜ.NH4
                NO3[:,:,:,t] = nutₜ.NO3
                PO4[:,:,:,t] = nutₜ.PO4
                DOC[:,:,:,t] = nutₜ.DOC
                DON[:,:,:,t] = nutₜ.DON
                DOP[:,:,:,t] = nutₜ.DOP
                POC[:,:,:,t] = nutₜ.POC
                PON[:,:,:,t] = nutₜ.PON
                POP[:,:,:,t] = nutₜ.POP
            end
        else
            nothing
        end
        model.nutrients = nutₜ;
        model.t += 1
    end
    if RunOption.NutOutputChoice == false
        write_nut_nc_alltime(g, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, RunParam.nTime)
    else
        return nothing
    end
end
"""
    PA_advect!(individuals, ΔT, vel)
Individual advection using a simple scheme
Particle sinking included
"""
function PA_advect!(model, ΔT, vel_itp)
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
    PA_TimeStep!(model, RunParam, velᵇ)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PA_TimeStep!(model::Model_struct, ΔT, velᵇ::velocity, resultspath)
    model.t += 1
    phyts_b,counts_p,consume_p=phyt_update(model, ΔT)
    model.individuals.phytos = phyts_b
    if model.individuals.zoos ≠ nothing
        zoos_b,counts_z,consume_z =zoo_update(model, ΔT)
        consume_p = sum_consume(consume_p, consume_z)
        model.individuals.zoos = zoos_b
        counts_p.graze = counts_z.graze
    end
    agent_num = size(phyts_b,1)
    write_pop_dynamics(model.t, agent_num, counts_p, resultspath)
    nutₜ,gtr = nut_update(model, velᵇ, consume_p, ΔT)
    model.nutrients = nutₜ
end
function PA_TimeStep!(model::Model_struct, ΔT, velᵇ::velocity)
    model.t += 1
    phyts_b,counts_p,consume_p=phyt_update(model, ΔT)
    model.individuals.phytos = phyts_b
    if model.individuals.zoos ≠ nothing
        zoos_b,counts_z,consume_z =zoo_update(model, ΔT)
        consume_p = sum_consume(consume_p, consume_z)
        model.individuals.zoos = zoos_b
        counts_p.graze = counts_z.graze
    end
    nutₜ,gtr = nut_update(model, velᵇ, consume_p, ΔT)
    model.nutrients = nutₜ
end

"""
    PA_advectRK4!(individuals, ΔT, vel_field, grid)
Individual advection using a RK4 method
Used for 3D double grids
Particle sinking not included
"""
function PA_advectRK4!(model, ΔT, vel_itps)
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
