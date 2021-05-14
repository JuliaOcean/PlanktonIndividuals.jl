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

function plankton(N::Int64, arch::Architecture, sp::Int64, params::Dict, maxN)
    rawdata = StructArray(x   = zeros(maxN), y   = zeros(maxN), z   = zeros(maxN),
                          iS  = zeros(maxN), Sz  = zeros(maxN), Bm  = zeros(maxN), 
                          Cq  = zeros(maxN), Nq  = zeros(maxN), Pq  = zeros(maxN), 
                          chl = zeros(maxN), gen = zeros(maxN), age = zeros(maxN), 
                          ac  = zeros(maxN), idx = zeros(maxN),
                          graz= zeros(maxN), mort= zeros(maxN), dvid= zeros(maxN),
                          xi  = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN)) 
    data = replace_storage(array_type(arch), rawdata)

    proc = StructArray(PS   = zeros(maxN), VDOC = zeros(maxN), VNH4 = zeros(maxN), VNO3 = zeros(maxN),
                       VPO4 = zeros(maxN), ρchl = zeros(maxN), resp = zeros(maxN), BS   = zeros(maxN), 
                       exu  = zeros(maxN), graz = zeros(maxN), mort = zeros(maxN), dvid = zeros(maxN))
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

function individuals(params::Dict, arch::Architecture, Nsp, N, maxN)
    plank_data=[]
    if Nsp > 9
        throw(ArgumentError("INDIVIDUALS: species must ≤ 9!"))
    else
        for i in 1:Nsp
            plank = plankton(N, arch, i, params, maxN)
            push!(plank_data, plank)
        end
        plank_name = plank_names[1:Nsp]
        planks = NamedTuple{plank_name}(plank_data)
        return individuals(planks,(;))
    end
end

function gen_individuals!(plank, N::Int64, g::AbstractGrid, arch::Architecture; mask = nothing)
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

    if mask ≠ nothing
        if size(mask) == (g.Nx, g.Ny, g.Nz)
            find_inds!(plank.data, g, arch)
            mask_individuals!(plank.data, mask, N, arch)
        else
            throw(ArgumentError("nut_mask: grid mismatch, size(mask) must equal to (grid.Nx, grid.Ny, grid.Nz)."))
        end
    end
end

@kernel function mask_individuals_kernel!(plank, mask)
    i = @index(Global)
    xi = plank.xi[i] - 2
    yi = plank.yi[i] - 2
    zi = plank.zi[i] - 2
    plank.ac[i] = mask[xi, yi, zi] * plank.ac[i]
end
function mask_individuals!(plank, mask, N, arch)
    kernel! = mask_individuals_kernel!(device(arch), 256, (N,))
    event = kernel!(plank, mask)
    wait(device(arch), event)
    return nothing
end