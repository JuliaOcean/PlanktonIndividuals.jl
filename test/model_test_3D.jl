using PlanktonIndividuals

grid = RegularRectilinearGrid(size = (16, 16, 16), spacing = (2, 2, 2), halo = (2, 2, 2))

model = PlanktonModel(CPU(), grid) 

TP = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.V)
TP = TP + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)

uvel = zeros(16,16,16,11)
vvel = zeros(16,16,16,11)
wvel = zeros(16,16,17,11)

for i in 1:11
    uvel[:,:,:,i] .= randn(16,16,16) .* 1e-4
    vvel[:,:,:,i] .= randn(16,16,16) .* 1e-4
    wvel[:,:,:,i] .= randn(16,16,17) .* 1e-4
end

# add boundary conditions for DOC
model.nutrients.DOC.bc.x.left  = 1.0e-3 # west boundary condition of 0.1 mmol/m^2/second
model.nutrients.DOC.bc.x.right = randn(20,20) .* 1e-3 # east boundary condition of -0.1 mmol/m^2/second
model.nutrients.DOC.bc.y.left  = randn(20,20,10) .* 1e-3 # south boundary condition of 0.1 mmol/m^2/second
model.nutrients.DOC.bc.y.right = randn(20,20,10) .* 1e-3 # north boundary condition of -0.1 mmol/m^2/second

sim = PlanktonSimulation(model, ΔT = 60, nΔT = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TPt = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.V)
TPt = TPt + sum(model.individuals.phytos.sp1.data.Pq .+ 
                model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)