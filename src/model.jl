function Model(grid, RunParam;
               individuals = setup_agents(RunParam,1.0,0.25,grid),
                 nutrients, 
                       PAR = read_default_IR_input(trunc(Int,RunParam.nTime*RunParam.DelT/3600)), 
                      temp = read_default_temp_input(trunc(Int,RunParam.nTime*RunParam.DelT/3600)), 
                    params = param_default,
                    output = create_output(individuals)
                   )
    return Model_struct(individuals, nutrients, grid, PAR, temp, params, output)
end

function ModelRun(model::Model_struct, RunParam::RunParams, RunOption::RunOptions)
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

    for t in 1:RunParam.nTime
        phyts_a = copy(model.individuals[t]) # read data from last time step
        phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, RunParam.DelT, phyts_a, model::Model_struct)
        if RunOption.VelChoice == false
            velᵇ = read_offline_vels(RunOption.VelOfflineOpt,trunc(Int,t*RunParam.DelT/3600))
        else
            velᵇ=store_vel[t]
        end
        if model.grid.Nx > 1 & model.grid.Ny > 1
            velᵈ = double_grid(velᵇ,model.grid)
            agent_move(phyts_b,velᵈ,model.grid,model.params["K_sink"],RunParam.DelT)
        elseif model.grid.Nx == 1 & model.grid.Ny == 1 & model.grid.Nz > 1
            agent_move(phyts_b,velᵇ,model.grid,model.params["K_sink"],RunParam.DelT) # for 1D only, use big grid velocities
        elseif model.grid.Nx == 1 & model.grid.Ny == 1 & model.grid.Nz == 1
            nothing #for 0D only
        end
        push!(model.individuals,phyts_b)
        write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,model.output)
        agent_num = size(phyts_b,1)
        nutₜ,gtr = nut_update(model, velᵇ, consume, RunParam.DelT, RunOption.Dim)
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
    end
    if RunOption.NutOutputChoice == false
        write_nut_nc_alltime(g, DIC, DIN, DOC, DON, POC, PON, RunParam.nTime)
    else
        return nothing
    end
end
