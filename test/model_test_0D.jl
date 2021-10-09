using PlanktonIndividuals, Serialization

grid = RegularRectilinearGrid(size = (1, 1, 1), spacing = (32, 32, 32), halo = (2, 2, 2))


phyt_params = deserialize(dirname(pathof(PlanktonIndividuals))*"/../test/param5.bin")

model = PlanktonModel(CPU(), grid;
                      N_species = 5,
                      N_individual = 1024,
                      max_individuals = 1024*10,
                      phyt_params = update_phyt_params(phyt_params, QuotaMode()))

TP = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.Δx .* grid.Δy .* grid.Δz)
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


sim = PlanktonSimulation(model, ΔT = 60, iterations = 10)

update!(sim)

TPt = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.Δx .* grid.Δy .* grid.Δz)
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
