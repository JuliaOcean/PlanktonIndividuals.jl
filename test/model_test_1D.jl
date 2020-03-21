using PlanktonIndividuals, Serialization
samples=dirname(pathof(PlanktonIndividuals))*"/../samples/"
RunOption=RunOptions(true, Dict(), true, Dict())
PhytoOpt = PlankOpt(1000, 2, Int(1e0), [1.8e-11, 1.8e-10], 1.0, 0.25)
RunParam=RunParams(10, 600, PhytoOpt, false, nothing)
g = deserialize(samples*"grid1D.bin");
store_vels = deserialize(samples*"uvw1D.bin");

nut_init = [2.0, 0.05,0.05,0.01,20.0, 0.0, 0.0, 0.0, 0.0, 0.0];
model = PI_Model(g, RunParam; nutrients = setup_nutrients(g,nut_init));

TP = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
         .* g.V)
TP = TP + sum(model.individuals.phytos[11,:])
for i in 1:10
    vel = store_vels[1]
    vel_itp = generate_vel_itp(model.grid, vel)
    PI_advect!(model, RunParam.ΔT, vel_itp)
    PI_TimeStep!(model, RunParam.ΔT, vel)
end

TPt = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
          .* g.V)
TPt = TPt + sum(model.individuals.phytos[11,:])
