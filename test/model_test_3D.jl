using PlanktonIndividuals, Serialization
samples=dirname(pathof(PlanktonIndividuals))*"/../samples/"
RunOption=RunOptions(true, Dict(), true, Dict())
g = deserialize(samples*"grid.bin");
store_vels = deserialize(samples*"uvw.bin");

nut_init = [2.0, 0.05,0.05,0.01,20.0, 0.0, 0.0, 0.0, 0.0, 0.0];
model = PI_Model(g, RunParam; nutrients = setup_nutrients(g,nut_init));

TP = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
         .* g.V)
TP = TP + sum(model.individuals.phytos[8,:] + model.individuals.phytos[5,:]*model.params["R_PC"])
for i in 1:10
    t = model.t
    vel = store_vels[i]
    vel_itp = generate_vel_itp(model.grid, vel)
    PI_advect!(model, RunParam.ΔT, vel_itp)
    PI_TimeStep!(model, RunParam.ΔT, vel)
end

TPt = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
          .* g.V)
TPt = TPt + sum(model.individuals.phytos[8,:] + model.individuals.phytos[5,:]*model.params["R_PC"])
