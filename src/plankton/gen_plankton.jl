mutable struct plankton
    data::AbstractArray
    proc::AbstractArray
    sp::Int64
    p::NamedTuple
end

struct individuals
    phytos::NamedTuple
    zoos::NamedTuple
end

const λ = 6

function plankton(N::Int64, arch::Architecture, sp::Int64, params::Dict)
    rawdata = StructArray(x   = zeros(λ*N), y   = zeros(λ*N), z   = zeros(λ*N),
                          iS  = zeros(λ*N), Sz  = zeros(λ*N), Bm  = zeros(λ*N), 
                          Cq  = zeros(λ*N), Nq  = zeros(λ*N), Pq  = zeros(λ*N), 
                          chl = zeros(λ*N), gen = zeros(λ*N), age = zeros(λ*N), 
                          ac  = zeros(λ*N), idx = zeros(λ*N),
                          graz= zeros(λ*N), mort= zeros(λ*N), dvid= zeros(λ*N),
                          xi  = zeros(Int,λ*N), yi  = zeros(Int,λ*N), zi  = zeros(Int,λ*N)) 
    data = replace_storage(array_type(arch), rawdata)

    proc = StructArray(PS   = zeros(λ*N), VDOC = zeros(λ*N), VNH4 = zeros(λ*N), VNO3 = zeros(λ*N),
                       VPO4 = zeros(λ*N), ρchl = zeros(λ*N), resp = zeros(λ*N), BS   = zeros(λ*N), 
                       exu  = zeros(λ*N), graz = zeros(λ*N), mort = zeros(λ*N), dvid = zeros(λ*N))
    proc_d = replace_storage(array_type(arch), proc)

    pkeys = collect(keys(params))
    tmp = zeros(length(param_names))
    for i in 1:length(param_names)
        if length(findall(x->x==string(param_names[i]),pkeys))==0
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(tmp)
    return plankton(data, proc_d, sp, p)
end

const plank_names=(:sp1, :sp2, :sp3, :sp4, :sp5, :sp6, :sp7, :sp8, :sp9)
const param_names=(:Nsuper, :Cquota, :mean, :var, :Chl2Cint, :α, :Φ, :Tempref, :TempAe, :TempCoeff,
                   :PCmax, :PC_b, :VDOCmax, :VDOC_b, :VNO3max, :VNH4max, :VN_b, :VPO4max, :VP_b,
                   :KsatDOC, :KsatNH4, :KsatNO3, :KsatPO4, :Cqmax, :Cqmin, :Nqmax, :Nqmin, :Pqmax, :Pqmin,
                   :Chl2N, :R_NC, :R_PC, :k_mtb, :k_mtb_b, :respir_a, :respir_b,
                   :grz_P, :dvid_type, :dvid_P, :dvid_stp, :dvid_reg, :dvid_stp2, :dvid_reg2,
                   :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP, :mortFracC, :mortFracN, :mortFracP)

function individuals(params::Dict, arch::Architecture)
    plank_data=[]
    Nsp = params["Nsp"]
    N = params["Nind"]
    if Nsp > 9
        throw(ArgumentError("INDIVIDUALS: species must ≤ 9!"))
    else
        for i in 1:Nsp
            plank = plankton(N,arch,i,params)
            push!(plank_data,plank)
        end
        plank_name = plank_names[1:Nsp]
        planks = NamedTuple{plank_name}(plank_data)
        return individuals(planks,(;))
    end
end

"""
    gen_didividuals!(plank, grid, arch)
Set up a series of agents following a normal distribution (mean,var)
'Nindivi' is agent number for each species, 'sp' is number of species, 'Nsuper' is the number of cells one agent represents,
'Cquota' is the initial biomass for one cell
'(x,y,z)' of an individual is the actual location not grid indices
"""
function gen_individuals!(plank, N::Int64, g::Grids, arch::Architecture)
    mean = plank.p.mean
    var = plank.p.var
    Cquota = plank.p.Cquota
    Nsuper = plank.p.Nsuper
    cqmax = plank.p.Cqmax
    cqmin = plank.p.Cqmin
    nqmax = plank.p.Nqmax
    nqmin = plank.p.Nqmin
    pqmax = plank.p.Pqmax
    pqmin = plank.p.Pqmin
    Chl2Cint = plank.p.Chl2Cint

    plank.data.x[1:N]   .= rand(rng_type(arch), N) .* (g.xF[g.Nx+g.Hx+1] - g.xF[g.Hx+1])                   # x
    plank.data.y[1:N]   .= rand(rng_type(arch), N) .* (g.yF[g.Ny+g.Hy+1] - g.yF[g.Hy+1])                   # y
    plank.data.z[1:N]   .= rand(rng_type(arch), N) .* (g.zF[g.Nz+g.Hz+1] - g.zF[g.Hz+1]) .+ g.zF[g.Hz+1]   # z
    plank.data.iS[1:N]  .= max.(1.0, randn(rng_type(arch), N) .* var .+ mean)                              # init_size
    plank.data.Sz[1:N]  .= copy(plank.data.iS[1:N])                                                        # size
    plank.data.Bm[1:N]  .= Cquota .* plank.data.Sz[1:N] .* Nsuper                                          # Bm
    plank.data.Cq[1:N]  .= rand(rng_type(arch), N) .* (cqmax - cqmin)  .+ cqmin                            # Cq
    plank.data.Cq[1:N]  .= plank.data.Cq[1:N] .* plank.data.Bm[1:N]                                        # Cq
    plank.data.Nq[1:N]  .= rand(rng_type(arch), N) .* (nqmax - nqmin)  .+ nqmin                            # Nq
    plank.data.Nq[1:N]  .= plank.data.Nq[1:N] .* plank.data.Bm[1:N]                                        # Nq
    plank.data.Pq[1:N]  .= rand(rng_type(arch), N) .* (pqmax - pqmin)  .+ pqmin                            # Pq
    plank.data.Pq[1:N]  .= plank.data.Pq[1:N] .* plank.data.Bm[1:N]                                        # Pq
    plank.data.chl[1:N] .= plank.data.Bm[1:N] .* Chl2Cint                                                  # Chl
    plank.data.gen[1:N] .= 1.0                                                                             # generation
    plank.data.age[1:N] .= 1.0                                                                             # age
    plank.data.ac[1:N]  .= 1.0                                                                             # active

    # if RunParam.Zoo == false
    #     return individuals(phyts0,nothing)
    # else
    #     zoos0 = gen_zooplkt(params, grid)
    #     return individuals(phyts0,zoos0)
    # end
end

# """
#     gen_zooplkt(ZooOpt, grid)
# Set up zooplankton individuals according to 'RunParam'
# """
# function gen_zooplkt(params, grid)
#     Nsp = params["Z_Nsp"]
#     N = params["Z_Nind"]
#     mean = params["Z_mean"]
#     var = params["Z_var"]
#     Cquota = params["Z_Cquota"]
#     Nsuper = params["Z_Nsuper"]
#     zoos0 = zeros(Real,10,N*Nsp)
#     zoos0[:,1] .= rand(Uniform(grid.xF[2],grid.xF[end]), N*Nsp)  # x
#     zoos0[:,2] .= rand(Uniform(grid.yF[2],grid.yF[end]), N*Nsp)  # y
#     zoos0[:,3] .= rand(Uniform(grid.zF[2],grid.zF[end-1]), N*Nsp)# z
#     zoos0[:,4] .= max.(1.0, rand(Normal(mean,var), N*Nsp))       # size
#     for i in 1:Nsp
#         lower = Int(1+(i-1)*N)
#         upper = Int(N+(i-1)*N)
#         zoos0[lower:upper,5] .= Cquota[i]*Nsuper                 # Bm
#         zoos0[:,8] .= i                                          # species
#     end
#     zoos0[:,5] .= zoos0[:,5] .* zoos0[:,4]                       # Bm
#     zoos0[:,6] .= 0.0                                            # Nq
#     zoos0[:,7] .= 0.0                                            # Pq
#     zoos0[:,9] .= 1.0                                            # generation
#     zoos0[:,10] .= 1.0                                           # age
#     return zoos0
# end
