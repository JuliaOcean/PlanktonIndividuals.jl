using PlanktonIndividuals, Serialization

grid = RectilinearGrid(size = (16, 16, 1), x = (0,32), y = (0,32), z = (0,-32))

model = PlanktonModel(CPU(), grid) 

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

TP = tot_mass(model.nutrients.PO4.data, grid) +
     tot_mass(model.nutrients.DOP.data, grid) +
     tot_mass(model.nutrients.POP.data, grid)
TP = TP + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)

uvel = zeros(16,16,1,11)
vvel = zeros(16,16,1,11)
wvel = zeros(16,16,2,11)

for i in 1:11
    uvel[:,:,:,i] .= randn(16,16,1) .* 1e-4
    vvel[:,:,:,i] .= randn(16,16,1) .* 1e-4
    wvel[:,:,:,i] .= randn(16,16,2) .* 1e-4
end

sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TPt = tot_mass(model.nutrients.PO4.data, grid) +
      tot_mass(model.nutrients.DOP.data, grid) +
      tot_mass(model.nutrients.POP.data, grid)
TPt = TPt + sum(model.individuals.phytos.sp1.data.Pq .+ 
                model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)

@testset "PlanktonIndividuals 2D tests:" begin
    @test isapprox(TP,TPt; atol=1e1)
end 