using PlanktonIndividuals

grid = RectilinearGrid(size = (1, 1, 16), x = (0,32), y = (0,32), z = (0,-32), 
                              topology = (Bounded, Bounded, Bounded), halo = (2, 2, 2))

model = PlanktonModel(CPU(), grid, mode = MacroMolecularMode()) 

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

uvel = zeros(2,1,16,11)
vvel = zeros(1,2,16,11)
wvel = zeros(1,1,17,11)

for i in 1:11
    uvel[:,:,:,i] .= randn(2,1,16) .* 1e-4
    vvel[:,:,:,i] .= randn(1,2,16) .* 1e-4
    wvel[:,:,:,i] .= randn(1,1,17) .* 1e-4
end

sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TPt = tot_mass(model.tracers.PO4.data, grid) +
      tot_mass(model.tracers.DOP.data, grid) +
      tot_mass(model.tracers.POP.data, grid)
TPt = TPt + sum(model.individuals.phytos.sp1.data.PST .+ 
               model.individuals.phytos.sp1.data.DNA .* model.individuals.phytos.sp1.p.R_PC_DNA .+
               model.individuals.phytos.sp1.data.RNA .* model.individuals.phytos.sp1.p.R_PC_RNA)


@testset "PlanktonIndividuals 1D tests:" begin
    @test isapprox(TP,TPt; atol=2e1)
end 