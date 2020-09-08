mutable struct individuals
    phytos::AbstractArray
    zoos::Union{Nothing,Array}
end
"""
    gen_agents(RunParams,grid)
Set up a series of agents following a normal distribution (mean,var)
'Nindivi' is agent number for each species, 'sp' is number of species, 'Nsuper' is the number of cells one agent represents,
'Cquota' is the initial biomass for one cell
'(x,y,z)' of an individual is the actual location not grid indices
"""
function gen_agents(::CPUs,RunParam::RunParams,grid)
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
    phyts0[1,:] .= rand(Uniform(grid.xF[grid.Hx+1], grid.xF[grid.Nx+grid.Hx+1]), N*Nsp)               # x
    phyts0[2,:] .= rand(Uniform(grid.yF[grid.Hy+1], grid.yF[grid.Ny+grid.Hy+1]), N*Nsp)               # y
    phyts0[3,:] .= rand(Uniform(grid.zF[grid.Hz+1], grid.zF[grid.Nz+grid.Hz+1]), N*Nsp)               # z
    phyts0[5,:] .= max.(1.0, rand(Normal(mean,var), N*Nsp))                    # size
    phyts0[4,:] .= copy(phyts0[5,:])                                           # init_size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[6,lower:upper]  .= Cquota[i]*Nsuper                             # Bm
        phyts0[11,lower:upper] .= i                                            # species
    end
    phyts0[6,:] .= phyts0[6,:] .* phyts0[5,:]                                  # Bm
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[7, lower:upper] .= copy(phyts0[6,lower:upper]) .* rand(Uniform(cqmin[i],cqmax[i]),N) # Cq
        phyts0[8, lower:upper] .= copy(phyts0[6,lower:upper]) .* rand(Uniform(nqmin[i],nqmax[i]),N) # Nq
        phyts0[9, lower:upper] .= copy(phyts0[6,lower:upper]) .* rand(Uniform(pqmin[i],pqmax[i]),N) # Pq
        phyts0[10,lower:upper] .= copy(phyts0[6,lower:upper]) .* params["Chl2Cint"][i]              # Chl
    end
    phyts0[12,:] .= 1.0                                                        # generation
    phyts0[13,:] .= 1.0                                                        # age

    if RunParam.Zoo == false
        return individuals(phyts0,nothing)
    else
        zoos0 = gen_zooplkt(params, grid)
        return individuals(phyts0,zoos0)
    end
end
function gen_agents(::GPUs,RunParam::RunParams,grid)
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
    phyts0[1,:] .= CuArray(rand(Uniform(grid.xF[grid.Hx+1], grid.xF[grid.Hx+grid.Nx+1]), N*Nsp))               # x
    phyts0[2,:] .= CuArray(rand(Uniform(grid.yF[grid.Hy+1], grid.yF[grid.Hy+grid.Ny+1]), N*Nsp))               # y
    phyts0[3,:] .= CuArray(rand(Uniform(grid.zF[grid.Hz+1], grid.zF[grid.Hz+grid.Nz+1]), N*Nsp))             # z
    phyts0[5,:] .= CuArray(max.(1.0, rand(Normal(mean,var), N*Nsp)))                    # size
    phyts0[4,:] .= copy(phyts0[5,:])                                                    # init_size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[6,lower:upper]  .= Cquota[i]*Nsuper                             # Bm
        phyts0[11,lower:upper] .= i                                            # species
    end
    phyts0[6,:] .= phyts0[6,:] .* phyts0[5,:]                                  # Bm
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        phyts0[7, lower:upper] .= copy(phyts0[6,lower:upper]) .* CuArray(rand(Uniform(cqmin[i],cqmax[i]),N)) # Cq
        phyts0[8, lower:upper] .= copy(phyts0[6,lower:upper]) .* CuArray(rand(Uniform(nqmin[i],nqmax[i]),N)) # Nq
        phyts0[9, lower:upper] .= copy(phyts0[6,lower:upper]) .* CuArray(rand(Uniform(pqmin[i],pqmax[i]),N)) # Pq
        phyts0[10,lower:upper] .= copy(phyts0[6,lower:upper]) .* params["Chl2Cint"][i]                       # Chl
    end
    phyts0[12,:] .= 1.0                                                        # generation
    phyts0[13,:] .= 1.0                                                        # age

    return individuals(phyts0,nothing)

    # if RunParam.Zoo == false
    #     return individuals(phyts0,nothing)
    # else
    #     zoos0 = setup_zooplkt(params, grid)
    #     return individuals(phyts0,zoos0)
    # end
end

"""
    gen_zooplkt(ZooOpt, grid)
Set up zooplankton individuals according to 'RunParam'
"""
function gen_zooplkt(params, grid)
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
