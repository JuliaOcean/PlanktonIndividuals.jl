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
    cqmax = params["Cqmax"]
    cqmin = params["Cqmin"]
    nqmax = params["Nqmax"]
    nqmin = params["Nqmin"]
    pqmax = params["Pqmax"]
    pqmin = params["Pqmin"]
    phyts0 = zeros(Real,13,N*Nsp)
    phyts0[1,:] .= rand(Uniform(grid.xF[2],grid.xF[end]), N*Nsp)         # x
    phyts0[2,:] .= rand(Uniform(grid.yF[2],grid.yF[end]), N*Nsp)         # y
    phyts0[3,:] .= rand(Uniform(grid.zF[2],grid.zF[end-1]), N*Nsp)       # z
    phyts0[4,:] .= max.(1.0, rand(Normal(mean,var), N*Nsp))              # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[5,lower:upper]  .= Cquota[i]*Nsuper                       # Bm
        phyts0[10,lower:upper] .= i                                      # species
    end
    phyts0[5,:] .= phyts0[5,:] .* phyts0[4,:]                            # Bm
    phyts0[6,:] .= copy(phyts0[5,:]) .* rand(Uniform(cqmin,cqmax),N*Nsp) # Cq
    phyts0[7,:] .= copy(phyts0[5,:]) .* rand(Uniform(nqmin,nqmax),N*Nsp) # Nq
    phyts0[8,:] .= copy(phyts0[5,:]) .* rand(Uniform(pqmin,pqmax),N*Nsp) # Pq
    phyts0[9,:] .= copy(phyts0[5,:]) .* params["Chl2Cint"]               # Chl
    phyts0[11,:] .= 1.0                                                  # generation
    phyts0[12,:] .= 1.0                                                  # age
    phyts0[13,:] .= copy(phyts0[4,:])                                    # init_size

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
    zoos0 = zeros(Real,10,N*Nsp)
    zoos0[1,:] .= rand(Uniform(grid.xF[2],grid.xF[end]), N*Nsp)  # x
    zoos0[2,:] .= rand(Uniform(grid.yF[2],grid.yF[end]), N*Nsp)  # y
    zoos0[3,:] .= rand(Uniform(grid.zF[2],grid.zF[end-1]), N*Nsp)# z
    zoos0[4,:] .= max.(1.0, rand(Normal(mean,var), N*Nsp))       # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        zoos0[5,lower:upper] .= Cquota[i]*Nsuper                 # Bm
        zoos0[8,:] .= i                                          # species
    end
    zoos0[5,:] .= zoos0[5,:] .* zoos0[4,:]                       # Bm
    zoos0[6,:] .= 0.0                                            # Nq
    zoos0[7,:] .= 0.0                                            # Pq
    zoos0[9,:] .= 1.0                                            # generation
    zoos0[10,:] .= 1.0                                           # age
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
'nut' is an array of 10 elements, each element is a kind of nutrient
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

"""
    load_nut_initials(paths,g)
Load nutrient initial conditions from files
"""
function load_nut_initials(paths,g)
    indices = ["DIC", "NH4", "NO3", "PO4", "DOC", "DON", "DOP", "POC", "PON", "POP"]
    pathkeys = collect(keys(paths))
    tmps = []
    for index in indices
        if length(findall(x->x==index, pathkeys)) == 0
            print("NUT_INIT: nutrient not found \n")
        else
            tmp = deserialize(paths[index])
            if size(tmp) == (g.Nx, g.Ny, g.Nz)
                push!(tmps,tmp)
            else
                print("NUT_INIT: grid mismatch \n")
            end
        end
    end
    nut = nutrient_fields(tmps[1],tmps[2],tmps[3],tmps[4],tmps[5],tmps[6],tmps[7],tmps[8],tmps[9],tmps[10])
    return nut
end
