function PA_Model(grid, RunParam;
                         t = 1,
               individuals = setup_agents(RunParam,1.0,0.25,grid),
                 nutrients, 
                       PAR = read_default_IR_input(RunParam.nTime, RunParam.DelT, grid), 
                      temp = read_default_temp_input(RunParam.nTime, RunParam.DelT, grid), 
                    params = param_default,
                    output = create_output(individuals)
                   )
    return Model_struct(t,individuals, nutrients, grid, PAR, temp, params, output)
end

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
        phyts_a = copy(model.individuals[t]) # read data from last time step
        phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, RunParam.DelT, phyts_a, model)
        if RunOption.VelChoice == false
            velᵇ = read_offline_vels(RunOption.VelOfflineOpt,trunc(Int,t*RunParam.DelT/3600))
        else
            velᵇ=store_vel[t]
        end
        if (model.grid.Nx > 1) & (model.grid.Ny > 1)
            velᵈ = double_grid_2D(velᵇ,model.grid)
            agent_advection(phyts_b,velᵈ,model.grid,model.params["k_sink"],RunParam.DelT)
        elseif (model.grid.Nx == 1) & (model.grid.Ny == 1) & (model.grid.Nz > 1)
            agent_advection(phyts_b,velᵇ,model.grid,model.params["k_sink"],RunParam.DelT) # for 1D only, use big grid velocities
        elseif (model.grid.Nx == 1) & (model.grid.Ny == 1) & (model.grid.Nz == 1)
            nothing #for 0D only
        end
        push!(model.individuals,phyts_b)
        write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,model.output)
        agent_num = size(phyts_b,1)
        nutₜ,gtr = nut_update(model, velᵇ, consume, RunParam.DelT, Dim = RunOption.Dim)
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

function PA_advect!(model::Model_struct, DelT, velᵇ::velocity)
    phyts_b = model.individuals[end]
    if (model.grid.Nx > 1) & (model.grid.Ny > 1)
        velᵈ = double_grid_2D(velᵇ)
        agent_advection(phyts_b,velᵈ,model.grid,model.params["k_sink"],DelT,"2D")
    elseif (model.grid.Nx == 1) & (model.grid.Ny == 1) & (model.grid.Nz > 1)
        agent_advection(phyts_b,velᵇ,model.grid,model.params["k_sink"],DelT,"1D") # for 1D only, use big grid velocities
    end
end


function PA_TimeStep!(model::Model_struct, DelT, velᵇ::velocity)
    phyts_a = copy(model.individuals[end]) # read data from last time step
    phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(model.t, DelT, phyts_a, model)
    push!(model.individuals,phyts_b)
    write_output(model.t,phyts_b,dvid_ct,graz_ct,death_ct,model.output)
    agent_num = size(phyts_b,1)
    nutₜ,gtr = nut_update(model, velᵇ, consume, DelT)
    model.nutrients = nutₜ
    model.t += 1
end

function PA_advectRK4!(model::Model_struct, DelT, vel_field)
    phyts_b = model.individuals[end] # read data from last time step
    if (model.grid.Nx > 1) & (model.grid.Ny > 1) & (model.grid.Nz > 1)
        vel_field_d = double_grid_3D.(vel_field)
        agent_advectionRK4(phyts_b,vel_field_d,model.grid,DelT,"3D")
    elseif (model.grid.Nx > 1) & (model.grid.Ny > 1) & (model.grid.Nz = 1)
        vel_field_d = double_grid_2D.(vel_field)
        agent_advectionRK4(phyts_b,vel_field_d,model.grid,DelT,"2D")
    elseif (model.grid.Nx == 1) & (model.grid.Ny == 1) & (model.grid.Nz > 1)
        agent_advectionRK4(phyts_b,vel_field,model.grid,DelT,"1D") # for 1D only, use big grid velocities
    end
end
