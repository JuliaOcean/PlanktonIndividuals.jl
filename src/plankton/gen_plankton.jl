#= index
1  2  3  4       5     6   7   8   9   10   11   12   13  14  15  16   17   18   19
x  y  z  size_i  size  Bm  Cq  Nq  Pq  chl  gen  age  xi  yi  zi  NH4  NO3  PO4  DOC

20  21        22  23    24    25    26    27    28    29  30   31   32    33    34  35  36
αI  TempFunc  PS  VDOC  VNH4  VNO3  VPO4  ρchl  resp  BS  exu  grz  mort  dvid  xt  yt  zt

37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57
u0  u1  v0  v1  w0  w1  xd  yd  zd  u1  v1  w1  u2  v2  w2  u3  v3  w3  u4  v4  w4
=#
struct plankton
    data::AbstractArray{Float64,2}
    sp::Int64
end

struct individuals
    phytos::NamedTuple
    zoos::NamedTuple
end

function plankton(N::Int64, arch::Architecture, sp)
    data  = zeros(Float64, N, 57) |> array_type(arch)
    return plankton(data, sp)
end

const plank_names=(:sp1, :sp2, :sp3, :sp4, :sp5, :sp6, :sp7, :sp8, :sp9)

function individuals(Nsp::Int64, N::Int64, arch::Architecture)
    plank_data = (plankton(N, arch, 1), plankton(N, arch, 2), plankton(N, arch, 3),
                  plankton(N, arch, 4), plankton(N, arch, 5), plankton(N, arch, 6),
                  plankton(N, arch, 7), plankton(N, arch, 8), plankton(N, arch, 9))
    if Nsp > 9
        throw(ArgumentError("INDIVIDUALS: species must ≤ 9!"))
    else
        plank_name = plank_names[1:Nsp]
        planks = NamedTuple{plank_name}(plank_data[1:Nsp])
        return individuals(planks,(;))
    end
end

"""
    gen_didividuals!(plank, RunParams, grid, arch)
Set up a series of agents following a normal distribution (mean,var)
'Nindivi' is agent number for each species, 'sp' is number of species, 'Nsuper' is the number of cells one agent represents,
'Cquota' is the initial biomass for one cell
'(x,y,z)' of an individual is the actual location not grid indices
"""
function gen_individuals!(plank, RunParam::RunParams, g::Grids, arch::Architecture)
    params = RunParam.params
    N = params["P_Nind"]
    mean = params["P_mean"][plank.sp]
    var = params["P_var"][plank.sp]
    Cquota = params["P_Cquota"][plank.sp]
    Nsuper = params["P_Nsuper"][plank.sp]
    cqmax = params["Cqmax"][plank.sp]
    cqmin = params["Cqmin"][plank.sp]
    nqmax = params["Nqmax"][plank.sp]
    nqmin = params["Nqmin"][plank.sp]
    pqmax = params["Pqmax"][plank.sp]
    pqmin = params["Pqmin"][plank.sp]

    plank.data[:,1] .= rand(Uniform(g.xF[g.Hx+1], g.xF[g.Nx+g.Hx+1]), N) |> array_type(arch)   # x
    plank.data[:,2] .= rand(Uniform(g.yF[g.Hy+1], g.yF[g.Ny+g.Hy+1]), N) |> array_type(arch)   # y
    plank.data[:,3] .= rand(Uniform(g.zF[g.Hz+1], g.zF[g.Nz+g.Hz+1]), N) |> array_type(arch)   # z
    plank.data[:,4] .= max.(1.0, rand(Normal(mean,var), N) |> array_type(arch))                # init_size
    plank.data[:,5] .= copy(plank.data[:,4])                                                   # size
    plank.data[:,6] .= Cquota .* plank.data[:,5]                                               # Bm
    plank.data[:,7] .= rand(Uniform(cqmin, cqmax), N) |> array_type(arch)                      # Cq
    plank.data[:,7] .= plank.data[:,7] .* plank.data[:,6]                                      # Cq
    plank.data[:,8] .= rand(Uniform(nqmin, nqmax), N) |> array_type(arch)                      # Nq
    plank.data[:,8] .= plank.data[:,8] .* plank.data[:,6]                                      # Nq
    plank.data[:,9] .= rand(Uniform(pqmin, pqmax), N) |> array_type(arch)                      # Pq
    plank.data[:,9] .= plank.data[:,9] .* plank.data[:,6]                                      # Pq
    plank.data[:,10] .= plank.data[:,6] .* params["Chl2Cint"][plank.sp]                        # Chl
    plank.data[:,11] .= 1.0                                                                    # generation
    plank.data[:,12] .= 1.0                                                                    # age

    # if RunParam.Zoo == false
    #     return individuals(phyts0,nothing)
    # else
    #     zoos0 = gen_zooplkt(params, grid)
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
    zoos0[:,1] .= rand(Uniform(grid.xF[2],grid.xF[end]), N*Nsp)  # x
    zoos0[:,2] .= rand(Uniform(grid.yF[2],grid.yF[end]), N*Nsp)  # y
    zoos0[:,3] .= rand(Uniform(grid.zF[2],grid.zF[end-1]), N*Nsp)# z
    zoos0[:,4] .= max.(1.0, rand(Normal(mean,var), N*Nsp))       # size
    for i in 1:Nsp
        lower = Int(1+(i-1)*N)
        upper = Int(N+(i-1)*N)
        zoos0[lower:upper,5] .= Cquota[i]*Nsuper                 # Bm
        zoos0[:,8] .= i                                          # species
    end
    zoos0[:,5] .= zoos0[:,5] .* zoos0[:,4]                       # Bm
    zoos0[:,6] .= 0.0                                            # Nq
    zoos0[:,7] .= 0.0                                            # Pq
    zoos0[:,9] .= 1.0                                            # generation
    zoos0[:,10] .= 1.0                                           # age
    return zoos0
end
