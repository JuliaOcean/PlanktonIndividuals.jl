using PlanktonIndividuals

grid = RectilinearGrid(size = (1, 1, 1), x = (0,32), y = (0,32), z = (0,-32))

model = PlanktonModel(CPU(), grid;
                      mode = ProteinMode(),
                      phyto = phyto_setup(nothing, [1024], 1),
                      max_individuals = 1024*10)

function tot_mass(tracer, g)
    mass = zeros(g.Nx, g.Ny, g.Nz)
    for i in 1:g.Nx
        for j in 1:g.Ny
            for k in 1:g.Nz
                mass[i,j,k] = tracer[i+g.Hx, j+g.Hy, k+g.Hz] * PlanktonIndividuals.Grids.volume(i+g.Hx, j+g.Hy, k+g.Hz, g)
            end
        end
    end
    return sum(mass)
end

TP = tot_mass(model.tracers.PO4.data, grid) +
     tot_mass(model.tracers.DOP.data, grid) +
     tot_mass(model.tracers.POP.data, grid)
TP = TP + sum(model.individuals.phytos.sp1.data.PST .+ 
              model.individuals.phytos.sp1.data.DNA .* model.individuals.phytos.sp1.p.R_PC_DNA .+
              model.individuals.phytos.sp1.data.RNA .* model.individuals.phytos.sp1.p.R_PC_RNA)


sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10)

update!(sim)

TPt = tot_mass(model.tracers.PO4.data, grid) +
      tot_mass(model.tracers.DOP.data, grid) +
      tot_mass(model.tracers.POP.data, grid)
TPt=TPt + sum(model.individuals.phytos.sp1.data.PST .+ 
              model.individuals.phytos.sp1.data.DNA .* model.individuals.phytos.sp1.p.R_PC_DNA .+
              model.individuals.phytos.sp1.data.RNA .* model.individuals.phytos.sp1.p.R_PC_RNA)

@testset "PlanktonIndividuals 0D ProteinMode tests:" begin
    @test isapprox(TP,TPt; atol=1e1)
end
