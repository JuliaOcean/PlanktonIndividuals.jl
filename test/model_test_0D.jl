using PlanktonIndividuals, Serialization
samples=dirname(pathof(PlanktonIndividuals))*"/../samples/"
RunOption=RunOptions(false, true, true, Dict(), true, Dict())
PhytoOpt = PlankOpt(1000, 2, Int(1e0), [1.8e-11, 1.8e-10], 1.0, 0.25)
RunParam=RunParams(10, 600, PhytoOpt, false, nothing)
g = deserialize(samples*"grid0D.bin");
nut_init = [2.0, 0.05,0.05,0.01,20.0, 0.0, 0.0, 0.0, 0.0, 0.0];
model = PI_Model(g, RunParam; nutrients = setup_nutrients(g,nut_init));

TP = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
         .* g.V)
TP = TP + sum(model.individuals.phytos[11,:])
for i in 1:10
    model.t += 1
    t = model.t
    phyts_b,counts_p,consume_p=PlanktonIndividuals.phyt_update(model, RunParam.ΔT)
    model.individuals.phytos = phyts_b
    nutₜ,gtr = PlanktonIndividuals.nut_update(model, consume_p, RunParam.ΔT)
    model.nutrients = nutₜ
end

TPt = sum((model.nutrients.PO4 .+ model.nutrients.DOP .+ model.nutrients.POP)
          .* g.V)
TPt = TPt + sum(model.individuals.phytos[11,:])
