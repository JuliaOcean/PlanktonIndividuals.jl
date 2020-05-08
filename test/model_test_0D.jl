using PlanktonIndividuals, Serialization
samples=dirname(pathof(PlanktonIndividuals))*"/../samples/"
RunOption=RunOptions(true, Dict(), true, Dict())
g = deserialize(samples*"grid0D.bin");
nut_init = [2.0, 0.05,0.05,0.01,20.0, 0.0, 0.0, 0.0, 0.0, 0.0];
model = PI_Model(g, RunParam; nutrients = setup_nutrients(g,nut_init));

TP = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
         .* g.V)
TP = TP + sum(model.individuals.phytos[8,:] + model.individuals.phytos[5,:]*model.params["R_PC"])
for i in 1:10
    model.t = model.t+RunParam.ΔT
    phyts_b,counts_p,consume_p=PlanktonIndividuals.phyt_update(model, RunParam.ΔT)
    model.individuals.phytos = phyts_b
    nutₜ,gtr = PlanktonIndividuals.nut_update(model, consume_p, RunParam.ΔT)
    model.nutrients = nutₜ
end

TPt = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
          .* g.V)
TPt = TPt + sum(model.individuals.phytos[8,:] + model.individuals.phytos[5,:]*model.params["R_PC"])
