function construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN)
    rawdata = StructArray(x    = zeros(maxN), y    = zeros(maxN), z    = zeros(maxN),
                          xi   = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN),
                          CH   = zeros(maxN), NST  = zeros(maxN), PST  = zeros(maxN),
                          PRO  = zeros(maxN), DNA  = zeros(maxN), RNA  = zeros(maxN), 
                          Chl  = zeros(maxN), gen  = zeros(maxN), age  = zeros(maxN), 
                          ac   = zeros(maxN), idx  = zeros(maxN),
                          PS   = zeros(maxN), VDOC = zeros(maxN), VNH4 = zeros(maxN),
                          VNO3 = zeros(maxN), VPO4 = zeros(maxN), ρChl = zeros(maxN),
                          resp = zeros(maxN), S_PRO= zeros(maxN), S_DNA= zeros(maxN),
                          S_RNA= zeros(maxN), exu  = zeros(maxN),
                          graz = zeros(maxN), mort = zeros(maxN), dvid = zeros(maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :C_DNA, :mean, :var, :PRO2DNA, :RNA2DNA, :Chl2DNA, :α, :Φ, :T⁺, :Ea,
                 :PCmax, :VDOCmax, :VNO3max, :VNH4max, :VPO4max, :KsatDOC, :KsatNH4, :KsatNO3,
                 :KsatPO4, :CHmax, :CHmin, :NSTmax, :NSTmin, :PSTmax, :PSTmin,
                 :Chl2N, :R_NC_PRO, :R_NC_DNA, :R_NC_RNA, :R_PC_DNA, :R_PC_RNA, :respir_a,
                 :k_pro_a, :k_sat_pro, :k_rna_a, :k_sat_rna, :k_dna_a, :k_sat_dna, 
                 :grz_P, :dvid_type, :dvid_P, :dvid_stp, :dvid_reg, :dvid_stp2, :dvid_reg2,
                 :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP, :mortFracC, :mortFracN, :mortFracP, :ther_mort)

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
    return plankton(data, p)
end

function generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)
    mean = plank.p.mean
    var = plank.p.var
    C_DNA = plank.p.C_DNA
    Nsuper = plank.p.Nsuper
    CHmax = plank.p.CHmax
    CHmin = plank.p.CHmin
    NSTmax = plank.p.NSTmax
    NSTmin = plank.p.NSTmin
    PSTmax = plank.p.PSTmax
    PSTmin = plank.p.PSTmin
    RNA2DNA = plank.p.RNA2DNA
    PRO2DNA = plank.p.PRO2DNA
    Chl2DNA = plank.p.Chl2DNA

    plank.data.ac[1:N]  .= 1.0                                                                             # activity
    plank.data.gen[1:N] .= 1.0                                                                             # generation
    plank.data.age[1:N] .= 0.0                                                                             # age

    randn!(rng_type(arch), plank.data.DNA)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.CH)
    rand!(rng_type(arch), plank.data.NST)
    rand!(rng_type(arch), plank.data.PST)

    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                         # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                         # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                         # z, unit: grid spacing, starting from 0
    plank.data.DNA .= max.(1.0, plank.data.DNA .* var .+ mean) .* C_DNA .* Nsuper .* plank.data.ac   # DNA
    plank.data.RNA .= plank.data.DNA .* RNA2DNA                                                      # RNA
    plank.data.PRO .= plank.data.DNA .* PRO2DNA                                                      # PRO
    plank.data.CH  .=(plank.data.CH  .* (CHmax  - CHmin)  .+ CHmin)  .* plank.data.PRO               # Cq
    plank.data.NST .=(plank.data.NST .* (NSTmax - NSTmin) .+ NSTmin) .* plank.data.PRO               # Nq
    plank.data.PST .=(plank.data.PST .* (PSTmax - PSTmin) .+ NSTmin) .* plank.data.PRO               # Pq
    plank.data.Chl .= plank.data.DNA .* Chl2DNA                                                      # Chl

    mask_individuals!(plank.data, g, N, arch)
end

@inline function functional_C_biomass(PRO, DNA, RNA)
    C_func = PRO + DNA + RNA
    return C_func
end
@inline function functional_N_biomass(PRO, DNA, RNA, p)
    N_func = PRO * p.R_NC_RPO + DNA * p.R_NC_DNA + RNA * p.R_NC_RNA
    return N_func
end
@inline function functional_P_biomass(DNA, RNA, p)
    P_func = DNA * p.R_PC_DNA + RNA * p.R_PC_RNA
    return P_func
end
@inline function total_C_biomass(PRO, DNA, RNA, CH)
    C_tot = PRO + DNA + RNA + CH
    return C_tot
end
@inline function total_N_biomass(PRO, DNA, RNA, NST, p)
    N_tot = PRO * p.R_NC_RPO + DNA * p.R_NC_DNA + RNA * p.R_NC_RNA + NST
    return N_tot
end
@inline function total_P_biomass(DNA, RNA, PST, p)
    P_tot = DNA * p.R_PC_DNA + RNA * p.R_PC_RNA + PST
    return P_tot
end