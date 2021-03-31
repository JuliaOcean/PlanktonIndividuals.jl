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

function plankton(N::Int64, arch::Architecture, sp::Int64, params::Dict, cap)
    rawdata = StructArray(x   = zeros(cap*N), y   = zeros(cap*N), z   = zeros(cap*N),
                          iS  = zeros(cap*N), Sz  = zeros(cap*N), Bm  = zeros(cap*N), 
                          Cq  = zeros(cap*N), Nq  = zeros(cap*N), Pq  = zeros(cap*N), 
                          chl = zeros(cap*N), gen = zeros(cap*N), age = zeros(cap*N), 
                          ac  = zeros(cap*N), idx = zeros(cap*N),
                          graz= zeros(cap*N), mort= zeros(cap*N), dvid= zeros(cap*N),
                          xi  = zeros(Int,cap*N), yi  = zeros(Int,cap*N), zi  = zeros(Int,cap*N)) 
    data = replace_storage(array_type(arch), rawdata)

    proc = StructArray(PS   = zeros(cap*N), VDOC = zeros(cap*N), VNH4 = zeros(cap*N), VNO3 = zeros(cap*N),
                       VPO4 = zeros(cap*N), ρchl = zeros(cap*N), resp = zeros(cap*N), BS   = zeros(cap*N), 
                       exu  = zeros(cap*N), graz = zeros(cap*N), mort = zeros(cap*N), dvid = zeros(cap*N))
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
const param_names=(:Nsuper, :Cquota, :mean, :var, :Chl2Cint, :α, :Φ, :TRef, :TAe, :TCoeff,
                   :PCmax, :PC_b, :VDOCmax, :VDOC_b, :VNO3max, :VNH4max, :VN_b, :VPO4max, :VP_b,
                   :KsatDOC, :KsatNH4, :KsatNO3, :KsatPO4, :Cqmax, :Cqmin, :Nqmax, :Nqmin, :Pqmax, :Pqmin,
                   :Chl2N, :R_NC, :R_PC, :k_mtb, :k_mtb_b, :respir_a, :respir_b,
                   :grz_P, :dvid_type, :dvid_P, :dvid_stp, :dvid_reg, :dvid_stp2, :dvid_reg2,
                   :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP, :mortFracC, :mortFracN, :mortFracP)

function individuals(params::Dict, arch::Architecture, Nsp, N, cap)
    plank_data=[]
    if Nsp > 9
        throw(ArgumentError("INDIVIDUALS: species must ≤ 9!"))
    else
        for i in 1:Nsp
            plank = plankton(N, arch, i, params, cap)
            push!(plank_data, plank)
        end
        plank_name = plank_names[1:Nsp]
        planks = NamedTuple{plank_name}(plank_data)
        return individuals(planks,(;))
    end
end

function gen_individuals!(plank, N::Int64, g::RegularRectilinearGrid, arch::Architecture)
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

    plank.data.ac[1:N]  .= 1.0                                                                             # activity
    plank.data.gen[1:N] .= 1.0                                                                             # generation
    plank.data.age[1:N] .= 0.0                                                                             # age

    randn!(rng_type(arch), plank.data.iS)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.Cq)
    rand!(rng_type(arch), plank.data.Nq)
    rand!(rng_type(arch), plank.data.Pq)

    plank.data.x   .=(plank.data.x .* (g.xF[g.Nx+g.Hx+1] - g.xF[g.Hx+1]) .+ g.xF[g.Hx+1]) .* plank.data.ac   # x
    plank.data.y   .=(plank.data.y .* (g.yF[g.Ny+g.Hy+1] - g.yF[g.Hy+1]) .+ g.yF[g.Hy+1]) .* plank.data.ac   # y
    plank.data.z   .=(plank.data.z .* (g.zF[g.Nz+g.Hz+1] - g.zF[g.Hz+1]) .+ g.zF[g.Hz+1]) .* plank.data.ac   # z
    plank.data.iS  .= max.(1.0, plank.data.iS .* var .+ mean) .* plank.data.ac                               # init_size
    plank.data.Sz  .= copy(plank.data.iS)                                                                    # size
    plank.data.Bm  .= Cquota .* plank.data.Sz .* Nsuper                                                      # Bm
    plank.data.Cq  .=(plank.data.Cq .* (cqmax - cqmin)  .+ cqmin) .* plank.data.Bm                           # Cq
    plank.data.Nq  .=(plank.data.Nq .* (nqmax - nqmin)  .+ nqmin) .* plank.data.Bm                           # Nq
    plank.data.Pq  .=(plank.data.Pq .* (pqmax - pqmin)  .+ pqmin) .* plank.data.Bm                           # Pq
    plank.data.chl .= plank.data.Bm .* Chl2Cint                                                              # Chl
end