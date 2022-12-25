using PlanktonIndividuals

grid = RectilinearGrid(size = (16, 16, 16), x = (0,32), y = (0,32), z = (0,-32))

model = PlanktonModel(CPU(), grid; mode = CarbonMode()) 

function tot_mass(nut, g)
    mass = zeros(g.Nx, g.Ny, g.Nz)
    for i in 1:g.Nx
        for j in 1:g.Ny
            for k in 1:g.Nz
                mass[i,j,k] = nut[i+g.Hx, j+g.Hy, k+g.Hz] * PlanktonIndividuals.Grids.volume(i+g.Hx, j+g.Hy, k+g.Hz, g)
            end
        end
    end
    return sum(mass)
end

TC = tot_mass(model.nutrients.DIC.data, grid) +
     tot_mass(model.nutrients.DOC.data, grid) +
     tot_mass(model.nutrients.POC.data, grid)
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
model.nutrients.DON.bc.east  = randn(16,16) .* 1e-3 # east boundary condition of -0.1 mmol/m^2/second
model.nutrients.DON.bc.south = randn(16,16,10) .* 1e-3 # south boundary condition of 0.1 mmol/m^2/second
model.nutrients.DON.bc.north = randn(16,16,10) .* 1e-3 # north boundary condition of -0.1 mmol/m^2/second

sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TCt = tot_mass(model.nutrients.DIC.data, grid) +
      tot_mass(model.nutrients.DOC.data, grid) +
      tot_mass(model.nutrients.POC.data, grid)
TCt = TCt + sum(model.individuals.phytos.sp1.data.Bm)

@testset "PlanktonIndividuals 3D tests:" begin
    @test isapprox(TC,TCt; atol=1e2)
end 