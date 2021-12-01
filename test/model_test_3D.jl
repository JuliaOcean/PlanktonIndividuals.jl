using PlanktonIndividuals

grid = RegularRectilinearGrid(size = (16, 16, 16), spacing = (2, 2, 2), halo = (2, 2, 2))

model = PlanktonModel(CPU(), grid; mode = CarbonMode()) 

TC = sum((interior(model.nutrients.DIC.data, grid) .+ 
          interior(model.nutrients.DOC.data, grid) .+ 
          interior(model.nutrients.POC.data, grid)) .* grid.Δx .* grid.Δy .* grid.Δz)
TC = TC + sum(model.individuals.phytos.sp1.data.Bm)

uvel = zeros(16,16,16,11)
vvel = zeros(16,16,16,11)
wvel = zeros(16,16,17,11)

for i in 1:11
    uvel[:,:,:,i] .= randn(16,16,16) .* 1e-4
    vvel[:,:,:,i] .= randn(16,16,16) .* 1e-4
    wvel[:,:,:,i] .= randn(16,16,17) .* 1e-4
end

# add boundary conditions for DOC
model.nutrients.DON.bc.west  = 1.0e-3 # west boundary condition of 0.1 mmol/m^2/second
model.nutrients.DON.bc.east  = randn(20,20) .* 1e-3 # east boundary condition of -0.1 mmol/m^2/second
model.nutrients.DON.bc.south = randn(20,20,10) .* 1e-3 # south boundary condition of 0.1 mmol/m^2/second
model.nutrients.DON.bc.north = randn(20,20,10) .* 1e-3 # north boundary condition of -0.1 mmol/m^2/second

sim = PlanktonSimulation(model, ΔT = 60, iterations = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TCt = sum((interior(model.nutrients.DIC.data, grid) .+ 
           interior(model.nutrients.DOC.data, grid) .+ 
           interior(model.nutrients.POC.data, grid)) .* grid.Δx .* grid.Δy .* grid.Δz)
TCt = TCt + sum(model.individuals.phytos.sp1.data.Bm)