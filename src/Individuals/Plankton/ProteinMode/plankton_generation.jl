function construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x      = zeros(FT, maxN), y      = zeros(FT, maxN), z      = zeros(FT, maxN),
                          xi     = zeros(Int,maxN), yi     = zeros(Int,maxN), zi     = zeros(Int,maxN),
                          CH     = zeros(FT, maxN), PST    = zeros(FT, maxN),
                          DNA    = zeros(FT, maxN), RNA    = zeros(FT, maxN), 
                          qFe    = zeros(FT, maxN), qNH4   = zeros(FT, maxN), qNO3   = zeros(FT, maxN),
                          PRO_RB = zeros(FT, maxN), PRO_MC = zeros(FT, maxN), PRO_MN = zeros(FT, maxN),
                          PRO_TN = zeros(FT, maxN), PRO_TP = zeros(FT, maxN), PRO_TFe= zeros(FT, maxN),
                          PRO_RS = zeros(FT, maxN), PRO_PS = zeros(FT, maxN),
                          Chl    = zeros(FT, maxN), gen    = zeros(FT, maxN), age    = zeros(FT, maxN), 
                          idx    = zeros(Int,maxN), ac     = zeros(Bool, maxN), 
                          PS     = zeros(FT, maxN), RS     = zeros(FT, maxN), ERS    = zeros(FT, maxN),
                          CF     = zeros(FT, maxN), NR     = zeros(FT, maxN), NF     = zeros(FT, maxN),
                          ECF    = zeros(FT, maxN), ENR    = zeros(FT, maxN), ENF    = zeros(FT, maxN),
                          EVNO3  = zeros(FT, maxN), EVPO4  = zeros(FT, maxN), EVFe   = zeros(FT, maxN),
                          ESP_RB = zeros(FT, maxN), ESP_MC = zeros(FT, maxN), ESP_MN = zeros(FT, maxN),
                          ESP_TN = zeros(FT, maxN), ESP_TP = zeros(FT, maxN), ESP_TFe= zeros(FT, maxN),
                          ESP_RS = zeros(FT, maxN), ESP_PS = zeros(FT, maxN), EDNA   = zeros(FT, maxN), 
                          ERNA   = zeros(FT, maxN), exE_RS = zeros(FT, maxN), exE_PS = zeros(FT, maxN),
                          VDOC   = zeros(FT, maxN), VNH4   = zeros(FT, maxN),
                          VNO3   = zeros(FT, maxN), VPO4   = zeros(FT, maxN), VFe    = zeros(FT, maxN),
                          ρChl   = zeros(FT, maxN), DChl   = zeros(FT, maxN),
                          DP_RB  = zeros(FT, maxN), DP_PS  = zeros(FT, maxN), DP_MC  = zeros(FT, maxN), 
                          SP_RB  = zeros(FT, maxN), SP_MC  = zeros(FT, maxN), SP_MN  = zeros(FT, maxN),
                          SP_TN  = zeros(FT, maxN), SP_TP  = zeros(FT, maxN), SP_TFe = zeros(FT, maxN),
                          SP_RS  = zeros(FT, maxN), SP_PS  = zeros(FT, maxN),
                          SDNA   = zeros(FT, maxN), SRNA   = zeros(FT, maxN), 
                          exu    = zeros(FT, maxN), ptc    = zeros(FT, maxN), Rptc   = zeros(FT, maxN),
                          graz   = zeros(FT, maxN), mort   = zeros(FT, maxN), dvid   = zeros(FT, maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :C_DNA, :var, :CH2DNA, :Chl2DNA,:RNA2DNA, 
                 :PRO_RB2DNA, :PRO_MC2DNA, :PRO_MN2DNA, :PRO_TN2DNA, 
                 :PRO_TP2DNA, :PRO_TFe2DNA, :PRO_RS2DNA, :PRO_PS2DNA, 
                 :α, :Topt, :Tmax, :Ea, :is_nr, :is_croc, :is_tric, :PCmax, :VDOCmax, 
                 :KcatCF, :KcatNF, :KcatNR, :KcatTNH4, :KcatTNO3, :KcatTPO4, :KcatTFe, :KcatRS, 
                 :r_max, :n_MC, :n_MN, :n_RB, :n_TN, :n_TPO4, :n_TFe, :n_RS, :n_PS,
                 :KsatDOC, :KsatNH4, :KsatNO3, :KsatPO4, :KsatNR, :KsatFe, :KFe_em,
                 :TN_max, :TP_max, :TFe_max,
                 :CHmax, :PSTmax, :qNH4max, :qNO3max, :qFemax,
                 :k_degRB, :k_degPS, :k_degMC, :k_degChl,
                 :PRO_RBmin, :PRO_PSmin, :PRO_MCmin, :Chlmin,
                 :Chl2N, :R_NC_PRO, :R_NC_DNA, :R_NC_RNA, :R_PC_DNA, :R_PC_RNA,:R_C_RNAPRB, 
                 :e_TNO3, :e_TPO4, :e_TFe, :e_RS, :e_CF, :e_NF, :e_NR, :e_SP, :e_DNA, :e_RNA, :e_min,
                 :k_sat_RNA, :k_DNA, :k_sat_DNA, :k_sat_PRO, 
                 :dvid_P, :grz_P, :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP, :grazFracFe,
                 :mortFracC, :mortFracN, :mortFracP, :mortFracFe)

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
    PRO_RB2DNA = plank.p.PRO_RB2DNA
    PRO_MC2DNA = plank.p.PRO_MC2DNA
    PRO_MN2DNA = plank.p.PRO_MN2DNA
    PRO_TN2DNA = plank.p.PRO_TN2DNA
    PRO_TP2DNA = plank.p.PRO_TP2DNA
    PRO_TFe2DNA= plank.p.PRO_TFe2DNA
    PRO_RS2DNA = plank.p.PRO_RS2DNA
    PRO_PS2DNA = plank.p.PRO_PS2DNA
    CH2DNA  = plank.p.CH2DNA
    Chl2DNA = plank.p.Chl2DNA

    plank.data.ac[1:N]  .= true                     # activity
    plank.data.gen[1:N] .= 1.0f0                    # generation
    plank.data.age[1:N] .= 0.0f0                    # age

    rand!(rng_type(arch), plank.data.DNA)
    rand!(rng_type(arch), plank.data.RNA)
    rand!(rng_type(arch), plank.data.PRO_RB)
    rand!(rng_type(arch), plank.data.PRO_MC)
    rand!(rng_type(arch), plank.data.PRO_MN)
    rand!(rng_type(arch), plank.data.PRO_TN)
    rand!(rng_type(arch), plank.data.PRO_TP)
    rand!(rng_type(arch), plank.data.PRO_TFe)
    rand!(rng_type(arch), plank.data.PRO_RS)
    rand!(rng_type(arch), plank.data.PRO_PS)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.CH)
    rand!(rng_type(arch), plank.data.PST)
    rand!(rng_type(arch), plank.data.qNO3)
    rand!(rng_type(arch), plank.data.qNH4)

    plank.data.DNA    .= plank.data.DNA    .* var .+ 1.0f0                  # range: (1.0,1.0+var)
    plank.data.RNA    .= plank.data.RNA    .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_RB .= plank.data.PRO_RB .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_MC .= plank.data.PRO_MC .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_MN .= plank.data.PRO_MN .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_TN .= plank.data.PRO_TN .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_TP .= plank.data.PRO_TP .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_TFe.= plank.data.PRO_TFe.* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_RS .= plank.data.PRO_RS .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO_PS .= plank.data.PRO_PS .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.CH     .= plank.data.CH     .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PST    .= plank.data.PST    .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.qNO3   .= plank.data.qNO3   .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.qNH4   .= plank.data.qNH4   .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)


    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                             # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                             # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                             # z, unit: grid spacing, starting from 0

    plank.data.DNA .= plank.data.DNA .* C_DNA .* Nsuper .* plank.data.ac                                 # DNA mmolC/individual
    plank.data.RNA .= plank.data.RNA .* C_DNA .* Nsuper .* plank.data.ac .* RNA2DNA                      # RNA mmolC/individual

    plank.data.PRO_RB .= plank.data.PRO_RB .* C_DNA .* Nsuper .* plank.data.ac .* PRO_RB2DNA
    plank.data.PRO_MC .= plank.data.PRO_MC .* C_DNA .* Nsuper .* plank.data.ac .* PRO_MC2DNA
    plank.data.PRO_MN .= plank.data.PRO_MN .* C_DNA .* Nsuper .* plank.data.ac .* PRO_MN2DNA
    plank.data.PRO_TN .= plank.data.PRO_TN .* C_DNA .* Nsuper .* plank.data.ac .* PRO_TN2DNA
    plank.data.PRO_TP .= plank.data.PRO_TP .* C_DNA .* Nsuper .* plank.data.ac .* PRO_TP2DNA
    plank.data.PRO_TFe.= plank.data.PRO_TFe.* C_DNA .* Nsuper .* plank.data.ac .* PRO_TFe2DNA
    plank.data.PRO_RS .= plank.data.PRO_RS .* C_DNA .* Nsuper .* plank.data.ac .* PRO_RS2DNA
    plank.data.PRO_PS .= plank.data.PRO_PS .* C_DNA .* Nsuper .* plank.data.ac .* PRO_PS2DNA
    plank.data.CH     .= plank.data.CH  .* plank.data.DNA .* CH2DNA                                         # CH  mmolC/individual
    plank.data.PST    .= plank.data.PST .* plank.data.CH ./ 106.0f0                                         # PST mmolP/individual
    plank.data.qNO3   .= plank.data.qNO3.* plank.data.CH ./ 106.0f0 .* 16.0f0                               # PST mmolP/individual
    plank.data.qNH4   .= plank.data.qNH4.* plank.data.CH ./ 106.0f0 .* 16.0f0                               # PST mmolP/individual
    plank.data.Chl    .= plank.data.DNA .* Chl2DNA * 893.49f0 / 55.0f0                                      # Chl mgChl/individual

    mask_individuals!(plank.data, g, N, arch)
end

@inline function total_C_biomass(PRO_RB, PRO_MC, PRO_MN, PRO_TN, PRO_TP, 
                                 PRO_TFe, PRO_RS, PRO_PS, DNA, RNA, CH, Chl)
    PRO = PRO_RB + PRO_MC + PRO_MN + PRO_TN + PRO_TP + PRO_TFe + PRO_RS + PRO_PS
    C_tot = PRO + DNA + RNA + CH + Chl / 893.49f0 * 55.0f0
    return C_tot
end
@inline function total_N_biomass(PRO_RB, PRO_MC, PRO_MN, PRO_TN, PRO_TP, 
                                 PRO_TFe, PRO_RS, PRO_PS, DNA, RNA, qNO3, qNH4, Chl, p)
    PRO = PRO_RB + PRO_MC + PRO_MN + PRO_TN + PRO_TP + PRO_TFe + PRO_RS + PRO_PS
    N_tot = PRO * p.R_NC_PRO + DNA * p.R_NC_DNA + RNA * p.R_NC_RNA + qNO3 + qNH4 + Chl / 893.49f0 * 4.0f0
    return N_tot
end
@inline function total_P_biomass(DNA, RNA, PST, p)
    P_tot = DNA * p.R_PC_DNA + RNA * p.R_PC_RNA + PST
    return P_tot
end
