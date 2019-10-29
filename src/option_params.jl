include(dirname(pathof(PhytoAgentModel))*"/param_default.jl")

# Options & Params
#                   Dim output, NutOutput, GridChoice, Gridoff, VelChoice, Veloff, SaveGrid, SaveVel, Test
RunOption=RunOptions(3, false,  true,      false,      Dict(),  false,     Dict(), false,    false,   false)
#                 nTime, DelT, Nindivi, Nsp, Nsuper,    Cquota(mmol/cell)
RunParam=RunParams(10,   3600, 100000,  2,   Int(1e15), [1.8e-11, 1.8e-10])

# Update parameter values based on .yaml file
function update_params(parameters::Dict, tmp::Dict)
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
