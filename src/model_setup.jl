"""
    setup_agents(RunParams,grid)
Set up a series of agents following a normal distribution (mean,var)
'Nindivi' is agent number for each species, 'sp' is number of species, 'Nsuper' is the number of cells one agent represents,
'Cquota' is the initial biomass for one cell
'(x,y,z)' of an individual is the actual location not grid indices
"""
function setup_agents(RunParam::RunParams,grid)
    params = RunParam.params
    Nsp = params["P_Nsp"]
    N = params["P_Nind"]
    mean = params["P_mean"]
    var = params["P_var"]
    Cquota = params["P_Cquota"]
    Nsuper = params["P_Nsuper"]
    phyts0 = zeros(Real,12,N*Nsp)
    phyts0[1,:]  = rand(Uniform(grid.xF[1],grid.xF[end]), N*Nsp)   # x
    phyts0[2,:]  = rand(Uniform(grid.yF[1],grid.yF[end]), N*Nsp)   # y
    phyts0[3,:]  = rand(Uniform(grid.zF[1],grid.zF[end]), N*Nsp)   # z
    phyts0[4,:] .= 1.0                                             # species
    phyts0[5,:] .= 1.0                                             # generation
    phyts0[6,:] .= 1.0                                             # age
    phyts0[7,:]  = max.(1.0, rand(Normal(mean,var), N*Nsp))        # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[8,lower:upper] .= Cquota[i]*Nsuper                  # Bm
    end
    phyts0[8,:]  = phyts0[8,:] .* phyts0[7,:]                      # Bm
    phyts0[9,:] .= 0.0                                             # Cq
    phyts0[10,:] .= 0.0                                            # Nq
    phyts0[11,:] .= 0.0                                            # Pq
    phyts0[12,:] = copy(phyts0[8,:]) .* 0.2 # Chl:C=0.2 gChl/molC  # Chl

    if RunParam.Zoo == false
        return individuals(phyts0,nothing)
    else
        zoos0 = setup_zooplkt(params, grid)
        return individuals(phyts0,zoos0)
    end
end

"""
    setup_zooplkt(ZooOpt, grid)
Set up zooplankton individuals according to 'RunParam'
"""
function setup_zooplkt(params, grid)
    Nsp = params["Z_Nsp"]
    N = params["Z_Nind"]
    mean = params["Z_mean"]
    var = params["Z_var"]
    Cquota = params["Z_Cquota"]
    Nsuper = params["Z_Nsuper"]
    zoos0 = zeros(Real,12,N*Nsp) # Cq = 0.0, chl = 0.0
    zoos0[1,:]  = rand(Uniform(grid.xF[1],grid.xF[end]), N*Nsp)  # x
    zoos0[2,:]  = rand(Uniform(grid.yF[1],grid.yF[end]), N*Nsp)  # y
    zoos0[3,:]  = rand(Uniform(grid.zF[1],grid.zF[end]), N*Nsp)  # z
    zoos0[4,:] .= 1.0                                            # species
    zoos0[5,:] .= 1.0                                            # generation
    zoos0[6,:] .= 1.0                                            # age
    zoos0[7,:]  = max.(1.0, rand(Normal(mean,var), N*Nsp))       # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        zoos0[8,lower:upper] .= Cquota[i]*Nsuper                 # Bm
    end
    zoos0[8,:]  = zoos0[8,:] .* zoos0[7,:]                       # Bm
    zoos0[10,:] = 0.0                                            # Nq
    zoos0[11,:] = 0.0                                            # Pq
    return zoos0
end


"""
    nutrients_init(g)
"""
function nutrients_init(g)
    nut = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    return nut
end


"""
    setup_nutrients(g,nut)
Set up initial nutrient fields according to grid information
'Nut' is an array of 10 elements, each element is a kind of nutrient
"""
function setup_nutrients(g,nut)
    DIC = fill(nut[1],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    NH4 = fill(nut[2],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    NO3 = fill(nut[3],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    PO4 = fill(nut[4],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    DOC = fill(nut[5],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    DON = fill(nut[6],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    DOP = fill(nut[7],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    POC = fill(nut[8],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    PON = fill(nut[9],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)
    POP = fill(nut[10],(g.Nx, g.Ny, g.Nz)) .* rand(Uniform(0.8,1.2), g.Nx, g.Ny, g.Nz)

    nutrients = nutrient_fields(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP)
    return nutrients
end

function pop_counts()
    return pop_counts(0,0,0)
end
