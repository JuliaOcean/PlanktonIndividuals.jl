using Oceananigans
using PlanktonIndividuals
using YAML
using Serialization, KernelDensity

# use Oceananigans to generate grid information
Nx = 15       # Number of grid points in x and y
Nz = 60       # Number of grid points in z
Δz = 2.0      # Grid spacing in x, y, z (meters)
grid = RegularCartesianGrid(size=(1,1,1), length=(Δz*Nx, Δz*Nx, Δz*Nz))
phy_grid = read_Ogrids(grid);

# set up time step etc.
RunParam.nTime = 1440*10
RunParam.ΔT = 60

# generate model structure
#           DIC  NH4   NO3   PO4   DOC  DON  DOP  POC   PON  POP
nut_init = [2.0, 0.50, 1.20, 0.20, 8.0, 0.5, 0.1, 10.0, 0.5, 0.1];
phy_model = PI_Model(phy_grid, RunParam;
                     nutrients = setup_nutrients(phy_grid,nut_init));
resultpath = PrepRunDir()

# generate an empty array to store cell size density
size_dens = zeros(25,RunParam.nTime)
xInd = collect(0.1:0.1:2.5)

# run the model for nTime time steps
for i in 1:RunParam.nTime
    PI_TimeStep!(phy_model, RunParam.ΔT, resultpath)
    ksd =kde(phy_model.individuals.phytos[7,:])
    iksd = InterpKDE(ksd)
    size_dens[:,i] = pdf(iksd,xInd)
end

# save cell size density results to binary file
open(resultpath*"size_dens.bin", "w") do io
    serialize(io, size_dens)
end
