"""
    PI_model(grid, RunParam)
Generate the model structure for time step
Default distribution of individuals is Normal distribution with 1.0 as mean and 0.25 as SD
Default PAR and temp are from ../samples
"""
function PI_Model(grid, RunParam;
                  t = 0,
                  individuals = setup_agents(RunParam,grid),
                  nutrients,
                  PAR = read_IR_input(RunParam.ΔT, grid),
                  temp = read_temp_input(RunParam.ΔT, grid),
                  params = RunParam.params,
                  diags = diags_setup(RunParam.nTime, RunParam.ΔT, grid, RunParam.params["diag_freq"], RunParam.params["diag_inds"])
                  )
    model = Model_struct(t,individuals, nutrients, grid, PAR, temp, params, diags)
    if RunParam.Zoo == true
        model.params["Grz_P"] = 0
    else
        nothing
    end
    return model
end
