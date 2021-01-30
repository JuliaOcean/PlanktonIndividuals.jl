using PlanktonIndividuals, Serialization

grid = gen_Grid(size = (1, 1, 1), spacing = (32, 32, 32), halo = (2, 2, 2))


phyt_params = deserialize(dirname(pathof(PlanktonIndividuals))*"/../test/param5.bin")

model = PI_Model(CPUs(), grid;
                 individual_size = (Nsp = 5, N = 1024, cap = 10),
                 phyt_params = update_phyt_params(phyt_params))

TP = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.V)
TP = TP + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC) + 
          sum(model.individuals.phytos.sp2.data.Pq .+ 
              model.individuals.phytos.sp2.data.Bm .* model.individuals.phytos.sp2.p.R_PC) +
          sum(model.individuals.phytos.sp3.data.Pq .+ 
              model.individuals.phytos.sp3.data.Bm .* model.individuals.phytos.sp3.p.R_PC) +
          sum(model.individuals.phytos.sp4.data.Pq .+ 
              model.individuals.phytos.sp4.data.Bm .* model.individuals.phytos.sp4.p.R_PC) +
          sum(model.individuals.phytos.sp5.data.Pq .+ 
              model.individuals.phytos.sp5.data.Bm .* model.individuals.phytos.sp5.p.R_PC)


sim = PI_simulation(model, ΔT = 60, nΔT = 10, diag_freq = 3600)

update!(sim)

TPt = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.V)
TPt=TPt + sum(model.individuals.phytos.sp1.data.Pq .+ 
                model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC) +
          sum(model.individuals.phytos.sp2.data.Pq .+ 
              model.individuals.phytos.sp2.data.Bm .* model.individuals.phytos.sp2.p.R_PC) +
          sum(model.individuals.phytos.sp3.data.Pq .+ 
              model.individuals.phytos.sp3.data.Bm .* model.individuals.phytos.sp3.p.R_PC) +
          sum(model.individuals.phytos.sp4.data.Pq .+ 
              model.individuals.phytos.sp4.data.Bm .* model.individuals.phytos.sp4.p.R_PC) +
          sum(model.individuals.phytos.sp5.data.Pq .+ 
              model.individuals.phytos.sp5.data.Bm .* model.individuals.phytos.sp5.p.R_PC)
