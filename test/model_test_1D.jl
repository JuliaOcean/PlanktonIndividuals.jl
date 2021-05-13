using PlanktonIndividuals, Serialization

grid = RegularRectilinearGrid(size = (1, 1, 16), spacing = (32, 32, 2), halo = (2, 2, 2))

model = PlanktonModel(CPU(), grid) 

TP = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.V)
TP = TP + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)

uvel = zeros(1,1,16,11)
vvel = zeros(1,1,16,11)
wvel = zeros(1,1,17,11)

for i in 1:11
    uvel[:,:,:,i] .= randn(1,1,16) .* 1e-4
    vvel[:,:,:,i] .= randn(1,1,16) .* 1e-4
    wvel[:,:,:,i] .= randn(1,1,17) .* 1e-4
end

sim = PlanktonSimulation(model, ΔT = 60, nΔT = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TPt = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.V)
TPt = TPt + sum(model.individuals.phytos.sp1.data.Pq .+ 
                model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)

