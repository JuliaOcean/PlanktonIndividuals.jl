# Options & Params
#                    GridChoice, Gridoff, VelChoice, Veloff
RunOption=RunOptions(false,      Dict(),  false,     Dict())

#                   Nindivi, Nsp, Nsuper,    Cquota(mmol/cell),  mean, var
PhytoOpt = PlankOpt(1000,    2,   Int(1e0),  [1.8e-11, 1.8e-10], 1.0,  0.25)

#                 Nindivi, Nsp, Nsuper,    Cquota(mmol/cell), mean, var
ZooOpt = PlankOpt(100,     1,   Int(1e0),  [1.8e-9],          1.0,  0.25)

#                  nTime, ΔT,   PhytoOpt, Zoo,   ZooOpt
RunParam=RunParams(10,    3600, PhytoOpt, false, ZooOpt)


"""
    update_params!(parameters, tmp)
Update parameter values based on .yaml file
'parameters' is default parameter set
'tmp' is the parameters need to update
"""
function update_params!(parameters::Dict, tmp::Dict)
    tmp_keys = collect(keys(tmp))
    for key in tmp_keys
        if length(findall(x->x==key, collect(keys(parameters))))==0
            print("parameter not found")
        else
            parameters[key] = tmp[key]
        end
    end
    return parameters
end
