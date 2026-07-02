function construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x    = zeros(FT, maxN), y    = zeros(FT, maxN), z    = zeros(FT, maxN),
                          xi   = zeros(Int,maxN), yi   = zeros(Int,maxN), zi   = zeros(Int,maxN),
                          CH   = zeros(FT, maxN), NST  = zeros(FT, maxN), PST  = zeros(FT, maxN),
                          DNA  = zeros(FT, maxN), RNA  = zeros(FT, maxN), 
                          PRO_R = zeros(FT, maxN), PRO_Mc = zeros(FT, maxN),  PRO_Mn = zeros(FT, maxN),
                          PRO_Tn = zeros(FT, maxN), PRO_Tp = zeros(FT, maxN), PRO_Tfe = zeros(FT, maxN),
                          PRO_C = zeros(FT, maxN),
                          Chl  = zeros(FT, maxN), gen  = zeros(FT, maxN), age  = zeros(FT, maxN), 
                          idx  = zeros(Int,maxN), ac   = zeros(Bool, maxN), 
                          PS   = zeros(FT, maxN), resp = zeros(FT, maxN),
                          CF   = zeros(FT, maxN), NR   = zeros(FT, maxN), NF   = zeros(FT, maxN),
                          VDOC = zeros(FT, maxN), VNH4 = zeros(FT, maxN),
                          VNO3 = zeros(FT, maxN), VPO4 = zeros(FT, maxN),  VFe  = zeros(FT, maxN),
                          ρChl = zeros(FT, maxN), qFe  = zeros(FT, maxN),
                          S_Pr= zeros(FT, maxN), S_Pm= zeros(FT, maxN), S_Ptno= zeros(FT, maxN),
                          S_Ptnh= zeros(FT, maxN), S_Ptp= zeros(FT, maxN), S_Ptfe= zeros(FT, maxN),
                          S_DNA= zeros(FT, maxN),S_RNA= zeros(FT, maxN), 
                          exu  = zeros(FT, maxN), ptc  = zeros(FT, maxN),
                          Rptc = zeros(FT, maxN),
                          graz = zeros(FT, maxN), mort = zeros(FT, maxN), dvid = zeros(FT, maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :C_DNA, :var, :CH2DNA, :Chl2DNA,:RNA2DNA, 
                 :PRO_R2DNA, :PRO_Mc2DNA, :PRO_Mn2DNA, :PRO_TN2DNA, :PRO_Tp2DNA, :PRO_Tfe2DNA, :PRO_C2DNA, 
                 :α, :Φ, :Topt, :Tmax, :Ea,
                 :PCmax, :VDOCmax, 
                 :KcatCF, :KcatNF, :KcatNR, :KcatTNH4, :KcatTNO3, :KcatTPO4; :KcatTFe, :KcatC, :KSAFe,
                 :r_max, :n_mC, :n_mN, :n_r, :n_tN, :n_tPO4, :n_tFe, :n_c,
                 :KsatDOC, :KsatNH4, :KsatNO3, :KsatPO4, :KsatFe, 
                 :CHmax, :NSTmax, :PSTmax,
                 :Chl2N, :R_NC_PRO, :R_NC_DNA, :R_NC_RNA, :R_PC_DNA, :R_PC_RNA,:R_C_RNAPr, 
                 :respir, :k_rna, :k_sat_rna, :k_dna, :k_sat_dna, 
                 :dvid_P, :grz_P, :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP,
                 :mortFracC, :mortFracN, :mortFracP)

    pkeys = Symbol.(collect(keys(params)))
    tmp = zeros(length(param_names))
    for i in 1:length(param_names)
        if param_names[i] ∉ pkeys
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(FT.(tmp))
    return phytoplankton(data, p)
end

function initialize_plankton!(plank, N::Int, g::AbstractGrid, arch::Architecture)
    var = plank.p.var
    C_DNA = plank.p.C_DNA
    Nsuper = plank.p.Nsuper
    RNA2DNA = plank.p.RNA2DNA
    PRO_R2DNA = plank.p.PRO_R2DNA
    PRO_Mc2DNA = plank.p.PRO_Mc2DNA
    PRO_Mn2DNA = plank.p.PRO_Mn2DNA
    PRO_TN2DNA = plank.p.PRO_TN2DNA
    PRO_Tp2DNA = plank.p.PRO_Tp2DNA
    PRO_Tfe2DNA = plank.p.PRO_Tfe2DNA
    PRO_C2DNA = plank.p.PRO_C2DNA
    CH2DNA  = plank.p.CH2DNA
    Chl2DNA = plank.p.Chl2DNA

    plank.data.ac[1:N]  .= true                     # activity
    plank.data.gen[1:N] .= 1.0f0                    # generation
    plank.data.age[1:N] .= 0.0f0                    # age

    rand!(rng_type(arch), plank.data.DNA)
    rand!(rng_type(arch), plank.data.RNA)
    rand!(rng_type(arch), plank.data.PRO_R)
    rand!(rng_type(arch), plank.data.PRO_Mc)
    rand!(rng_type(arch), plank.data.PRO_Mn)
    rand!(rng_type(arch), plank.data.PRO_Tn)
    rand!(rng_type(arch), plank.data.PRO_Tp)
    rand!(rng_type(arch), plank.data.PRO_Tfe)
    rand!(rng_type(arch), plank.data.PRO_C)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.CH)
    rand!(rng_type(arch), plank.data.NST)
    rand!(rng_type(arch), plank.data.PST)

    plank.data.DNA .= plank.data.DNA .* var .+ 1.0f0                  # range: (1.0,1.0+var)
    plank.data.RNA .= plank.data.RNA .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_R .= plank.data.PRO_R .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_Mc .= plank.data.PRO_Mc .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_Mn .= plank.data.PRO_Mn .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_Tn .= plank.data.PRO_Tn .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_Tp .= plank.data.PRO_Tp .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_Tfe .= plank.data.PRO_Tfe .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_C .= plank.data.PRO_C .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.CH  .= plank.data.CH  .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.NST .= plank.data.NST .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PST .= plank.data.PST .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)


    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                             # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                             # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                             # z, unit: grid spacing, starting from 0
    plank.data.DNA .= plank.data.DNA .* C_DNA .* Nsuper .* plank.data.ac                                 # DNA mmolC/individual
    plank.data.RNA .= plank.data.RNA .* C_DNA .* Nsuper .* plank.data.ac .* RNA2DNA                      # RNA mmolC/individual
    plank.data.PRO_R .= plank.data.PRO_R .* C_DNA .* Nsuper .* plank.data.ac .* PRO_R2DNA
    plank.data.PRO_Mc .= plank.data.PRO_Mc .* C_DNA .* Nsuper .* plank.data.ac .* PRO_Mc2DNA
    plank.data.PRO_Mn .= plank.data.PRO_Mn .* C_DNA .* Nsuper .* plank.data.ac .* PRO_Mn2DNA
    plank.data.PRO_TN .= plank.data.PRO_TN .* C_DNA .* Nsuper .* plank.data.ac .* PRO_TN2DNA
    plank.data.PRO_Tp .= plank.data.PRO_Tp .* C_DNA .* Nsuper .* plank.data.ac .* PRO_Tp2DNA
    plank.data.PRO_Tfe .= plank.data.PRO_Tfe .* C_DNA .* Nsuper .* plank.data.ac .* PRO_Tfe2DNA
    plank.data.PRO_C .= plank.data.PRO_C .* C_DNA .* Nsuper .* plank.data.ac .* PRO_C2DNA
    plank.data.CH  .= plank.data.CH  .* plank.data.DNA .* CH2DNA                                         # CH  mmolC/individual
    plank.data.NST .= plank.data.NST .* plank.data.CH .* 16.0f0 ./ 106.0f0                               # NST mmolN/individual
    plank.data.PST .= plank.data.PST .* plank.data.CH ./ 106.0f0                                         # PST mmolP/individual
    plank.data.Chl .= plank.data.DNA .* Chl2DNA * 893.49f0 / 55.0f0                                      # Chl mgChl/individual

    mask_individuals!(plank.data, g, N, arch)
end

@inline function total_C_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_TN, PRO_Tp, PRO_Tfe, DNA, RNA, CH, Chl)
    PRO = PRO_R + PRO_Mc + PRO_Mn + PRO_TN + PRO_Tp + PRO_Tfe
    C_tot = PRO + DNA + RNA + CH + Chl / 893.49f0 * 55.0f0 
    return C_tot
end
@inline function total_N_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_TN, PRO_Tp, PRO_Tfe, DNA, RNA, NST, Chl, p)
    PRO = PRO_R + PRO_Mc + PRO_Mn + PRO_TN + PRO_Tp + PRO_Tfe
    N_tot = PRO * p.R_NC_PRO + DNA * p.R_NC_DNA + RNA * p.R_NC_RNA + NST + Chl / 893.49f0 * 4.0f0
    return N_tot
end
@inline function total_P_biomass(DNA, RNA, PST, p)
    P_tot = DNA * p.R_PC_DNA + RNA * p.R_PC_RNA + PST
    return P_tot
end
