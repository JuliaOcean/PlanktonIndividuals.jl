"""
    setup_agents(RunParams,grid)
Set up a series of agents following a normal distribution (mean,var)
'Nindivi' is agent number for each species, 'sp' is number of species, 'Nsuper' is the number of cells one agent represents,
'Cquota' is the initial biomass for one cell
'(x,y,z)' of an individual is the actual location not grid indices
"""
function setup_agents(::CPUs,RunParam::RunParams,grid)
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
    phyts0 = zeros(Float64,13,N*Nsp)
    phyts0[1,:] .= rand(Uniform(grid.xF[2],grid.xF[end]), N*Nsp)               # x
    phyts0[2,:] .= rand(Uniform(grid.yF[2],grid.yF[end]), N*Nsp)               # y
    phyts0[3,:] .= rand(Uniform(grid.zF[2],grid.zF[end-1]), N*Nsp)             # z
    phyts0[4,:] .= max.(1.0, rand(Normal(mean,var), N*Nsp))                    # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[5,lower:upper]  .= Cquota[i]*Nsuper                             # Bm
        phyts0[10,lower:upper] .= i                                            # species
    end
    phyts0[5,:] .= phyts0[5,:] .* phyts0[4,:]                                  # Bm
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[6,lower:upper] .= copy(phyts0[5,lower:upper]) .* rand(Uniform(cqmin[i],cqmax[i]),N) # Cq
        phyts0[7,lower:upper] .= copy(phyts0[5,lower:upper]) .* rand(Uniform(nqmin[i],nqmax[i]),N) # Nq
        phyts0[8,lower:upper] .= copy(phyts0[5,lower:upper]) .* rand(Uniform(pqmin[i],pqmax[i]),N) # Pq
        phyts0[9,lower:upper] .= copy(phyts0[5,lower:upper]) .* params["Chl2Cint"][i]              # Chl
    end
    phyts0[11,:] .= 1.0                                                        # generation
    phyts0[12,:] .= 1.0                                                        # age
    phyts0[13,:] .= copy(phyts0[4,:])                                          # init_size

    if RunParam.Zoo == false
        return individuals(phyts0,nothing)
    else
        zoos0 = setup_zooplkt(params, grid)
        return individuals(phyts0,zoos0)
    end
end
function setup_agents(::GPUs,RunParam::RunParams,grid)
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
    phyts0 = zeros(Float64,13,N*Nsp) |> CuArray
    phyts0[1,:] .= CuArray(rand(Uniform(grid.xF[2],grid.xF[end]), N*Nsp))               # x
    phyts0[2,:] .= CuArray(rand(Uniform(grid.yF[2],grid.yF[end]), N*Nsp))               # y
    phyts0[3,:] .= CuArray(rand(Uniform(grid.zF[2],grid.zF[end-1]), N*Nsp))             # z
    phyts0[4,:] .= CuArray(max.(1.0, rand(Normal(mean,var), N*Nsp)))                    # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[5,lower:upper]  .= Cquota[i]*Nsuper                             # Bm
        phyts0[10,lower:upper] .= i                                            # species
    end
    phyts0[5,:] .= phyts0[5,:] .* phyts0[4,:]                                  # Bm
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[6,lower:upper] .= copy(phyts0[5,lower:upper]) .* CuArray(rand(Uniform(cqmin[i],cqmax[i]),N)) # Cq
        phyts0[7,lower:upper] .= copy(phyts0[5,lower:upper]) .* CuArray(rand(Uniform(nqmin[i],nqmax[i]),N)) # Nq
        phyts0[8,lower:upper] .= copy(phyts0[5,lower:upper]) .* CuArray(rand(Uniform(pqmin[i],pqmax[i]),N)) # Pq
        phyts0[9,lower:upper] .= copy(phyts0[5,lower:upper]) .* params["Chl2Cint"][i]                       # Chl
    end
    phyts0[11,:] .= 1.0                                                        # generation
    phyts0[12,:] .= 1.0                                                        # age
    phyts0[13,:] .= copy(phyts0[4,:])                                          # init_size

    return individuals(phyts0,nothing)

    # if RunParam.Zoo == false
    #     return individuals(phyts0,nothing)
    # else
    #     zoos0 = setup_zooplkt(params, grid)
    #     return individuals(phyts0,zoos0)
    # end
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
function nutrients_init(::CPUs,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    nut = nutrient_fields(zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size))
    return nut
end
function nutrients_init(::GPUs,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    nut = nutrient_fields(CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)))
    return nut
end


"""
    setup_nutrients(g,nut)
Set up initial nutrient fields according to grid information
'nut' is an array of 10 elements, each element is a kind of nutrient
"""
function setup_nutrients(::CPUs,g,nut)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    DIC = fill(nut[1],total_size) .* rand(Uniform(0.8,1.2), total_size)
    NH4 = fill(nut[2],total_size) .* rand(Uniform(0.8,1.2), total_size)
    NO3 = fill(nut[3],total_size) .* rand(Uniform(0.8,1.2), total_size)
    PO4 = fill(nut[4],total_size) .* rand(Uniform(0.8,1.2), total_size)
    DOC = fill(nut[5],total_size) .* rand(Uniform(0.8,1.2), total_size)
    DON = fill(nut[6],total_size) .* rand(Uniform(0.8,1.2), total_size)
    DOP = fill(nut[7],total_size) .* rand(Uniform(0.8,1.2), total_size)
    POC = fill(nut[8],total_size) .* rand(Uniform(0.8,1.2), total_size)
    PON = fill(nut[9],total_size) .* rand(Uniform(0.8,1.2), total_size)
    POP = fill(nut[10],total_size) .* rand(Uniform(0.8,1.2), total_size)

    nutrients = nutrient_fields(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP)
    return nutrients
end
function setup_nutrients(::GPUs,g,nut)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    DIC = fill(nut[1],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    NH4 = fill(nut[2],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    NO3 = fill(nut[3],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    PO4 = fill(nut[4],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    DOC = fill(nut[5],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    DON = fill(nut[6],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    DOP = fill(nut[7],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    POC = fill(nut[8],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    PON = fill(nut[9],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    POP = fill(nut[10],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray

    nutrients = nutrient_fields(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP)
    return nutrients
end

"""
    load_nut_initials(paths,g)
Load nutrient initial conditions from files
"""
function load_nut_initials(::CPUs,paths,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    indices = ["DIC", "NH4", "NO3", "PO4", "DOC", "DON", "DOP", "POC", "PON", "POP"]
    pathkeys = collect(keys(paths))
    tmps = []
    for index in indices
        if length(findall(x->x==index, pathkeys)) == 0
            print("NUT_INIT: nutrient not found \n")
        else
            tmp = deserialize(paths[index])
            if size(tmp) == total_size
                push!(tmps,tmp)
            else
                print("NUT_INIT: grid mismatch \n")
            end
        end
    end
    nut = nutrient_fields(tmps[1],tmps[2],tmps[3],tmps[4],tmps[5],tmps[6],tmps[7],tmps[8],tmps[9],tmps[10])
    return nut
end
function load_nut_initials(::GPUs,paths,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    indices = ["DIC", "NH4", "NO3", "PO4", "DOC", "DON", "DOP", "POC", "PON", "POP"]
    pathkeys = collect(keys(paths))
    tmps = []
    for index in indices
        if length(findall(x->x==index, pathkeys)) == 0
            print("NUT_INIT: nutrient not found \n")
        else
            tmp = deserialize(paths[index]) |> CuArray
            if size(tmp) == total_size
                push!(tmps,tmp)
            else
                print("NUT_INIT: grid mismatch \n")
            end
        end
    end
    nut = nutrient_fields(tmps[1],tmps[2],tmps[3],tmps[4],tmps[5],tmps[6],tmps[7],tmps[8],tmps[9],tmps[10])
    return nut
end
