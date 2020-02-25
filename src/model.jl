"""
    PA_model(grid, RunParam)
Generate the model structure for time step
Default distribution of individuals is Normal distribution with 1.0 as mean and 0.25 as SD
Default PAR and temp are from ../samples
"""
function PA_Model(grid, RunParam;
                         t = 1,
               individuals = setup_agents(RunParam,grid),
                 nutrients,
                       PAR = read_default_IR_input(RunParam.nTime, RunParam.ΔT, grid),
                      temp = read_default_temp_input(RunParam.nTime, RunParam.ΔT, grid),
                    params = param_default,
                    output = create_output(individuals[:,1])
                   )
    return Model_struct(t,individuals, nutrients, grid, PAR, temp, params, output)
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
        DIN = zeros(g.Nx, g.Ny, g.Nz, nTime)
        DOC = zeros(g.Nx, g.Ny, g.Nz, nTime)
        DON = zeros(g.Nx, g.Ny, g.Nz, nTime)
        POC = zeros(g.Nx, g.Ny, g.Nz, nTime)
        PON = zeros(g.Nx, g.Ny, g.Nz, nTime)
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
            velᵇ=store_vel[model.t]
        end
        if (model.grid.Nx == 1) & (model.grid.Ny == 1) & (model.grid.Nz == 1)
            nothing #for 0D only
        else
            PA_advect!(model, RunParam.ΔT, velᵇ)
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
                DIN[:,:,:,t] = nutₜ.DIN
                DOC[:,:,:,t] = nutₜ.DOC
                DON[:,:,:,t] = nutₜ.DON
                POC[:,:,:,t] = nutₜ.POC
                PON[:,:,:,t] = nutₜ.PON
            end
        else
            nothing
        end
        model.nutrients = nutₜ;
        model.t += 1
    end
    if RunOption.NutOutputChoice == false
        write_nut_nc_alltime(g, DIC, DIN, DOC, DON, POC, PON, RunParam.nTime)
    else
        return nothing
    end
end
"""
    PA_advect!(individuals, ΔT, velᵇ)
Individual advection using a simple scheme
Used for 2D and 1D double grids
Particle sinking included
"""
function PA_advect!(model, ΔT, velᵇ::velocity)
    individuals = model.individuals
    grid = model.grid
    params = model.params
    for i in 1:size(individuals,2)
        if (grid.Nx > 1) & (grid.Ny > 1)
            velᵈ = double_grid_2D(velᵇ)
            agent_advection(individuals[end,i],velᵈ,grid,ΔT,"2D")
            agent_diffusionH(individuals[end,i],grid,params["κhP"])
            if grid.Nz > 1
                agent_diffusionV(individuals[end,i],grid,params["κvP"])
            else
                nothing
            end
        elseif (grid.Nx == 1) & (grid.Ny == 1) & (grid.Nz > 1)
            agent_advection(individuals[end,i],velᵇ,grid,ΔT,"1D") # for 1D only, use big grid velocities
            agent_diffusionV(individuals[end,i],grid,params["κvP"])
        end
    end
end

"""
    PA_TimeStep!(model, RunParam, velᵇ)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PA_TimeStep!(model::Model_struct, ΔT, velᵇ::velocity)
    phyts_a = copy(model.individuals[end,1]) # read data from last time step
    phyts_b,counts_p,consume_p=phyt_update(model.t, ΔT, phyts_a, model)
    if size(model.individuals,2) == 2
        zoos_a  = copy(model.individuals[end,2]) # read data from last time step
        zoos_b,counts_z,consume_z =zoo_update(zoos_a, phyts_b, ΔT, model)
        consume_p = sum_consume(consume_p, consume_z)
        individualsₜ = hcat([phyts_b],[zoos_b])
        model.individuals = vcat(model.individuals, individualsₜ)
    elseif size(model.individuals,2) == 1
        individualsₜ = [phyts_b]
        model.individuals = vcat(model.individuals, individualsₜ)
    end
    write_output(model.t,phyts_b,counts_p,model.output)
    agent_num = size(phyts_b,1)
    nutₜ,gtr = nut_update(model, velᵇ, consume_p, ΔT)
    model.nutrients = nutₜ
    model.t += 1
end

"""
    PA_advectRK4!(individuals, ΔT, vel_field, grid)
Individual advection using a RK4 method
Used for 3D double grids
Particle sinking not included
"""
function PA_advectRK4!(model, ΔT, vel_field)
    individuals = model.individuals
    grid = model.grid
    params = model.params
    for i in 1:size(individuals,2)
        if (grid.Nx > 1) & (grid.Ny > 1) & (grid.Nz > 1)
            vel_field_d = double_grid_3D.(vel_field)
            agent_advectionRK4(individuals[end,i],vel_field_d,grid,ΔT,"3D")
            agent_diffusionH(individuals[end,i],grid,params["κhP"])
            agent_diffusionV(individuals[end,i],grid,params["κvP"])
        elseif (grid.Nx > 1) & (grid.Ny > 1) & (grid.Nz == 1)
            vel_field_d = double_grid_2D.(vel_field)
            agent_advectionRK4(individuals[end,i],vel_field_d,grid,ΔT,"2D")
            agent_diffusionH(individuals[end,i],grid,params["κhP"])
        elseif (grid.Nx == 1) & (grid.Ny == 1) & (grid.Nz > 1)
            agent_advectionRK4(individuals[end,i],vel_field,grid,ΔT,"1D") # for 1D only, use big grid velocities
            agent_diffusionV(individuals[end,i],grid,params["κvP"])
        end
    end
end
