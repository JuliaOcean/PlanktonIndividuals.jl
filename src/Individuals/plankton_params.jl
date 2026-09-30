# Individual parameter defaults, dispatched by physiology mode.
function phyt_params_default end
function colony_params_default end
function abiotic_params_default end

function generate_n_species_params(N, params)
    p = []
    for key in keys(params)
        pa = (key, fill(params[key][1],N))
        push!(p, pa)
    end
    return Dict(p)
end

#=
CH includes cabohydrate and lipids
┌─────────┬────────────────┬────────────────────┬───────────────────────────┬──────────────────────────┐
│         │ Micromonas sp. │ Ostreococcus tauri │ Thalassiosira weissflogii │ Thalassiosira pseudonana │
├─────────┼────────────────┼────────────────────┼───────────────────────────┼──────────────────────────┤
│ PRO:DNA │          34.83 │              22.75 │                     67.99 │                   135.36 │
│ RNA:DNA │           1.73 │               1.71 │                      3.70 │                     6.49 │
│ CHL:DNA │           3.45 │               2.16 │                      9.83 │                    18.05 │
│  CH:DNA │          23.32 │              19.01 │                     79.26 │                   121.48 │
│   k_pro │       7.06e-05 │           1.50e-04 │                  7.64e-05 │                 1.17e-04 │
│   k_dna │       1.62e-07 │           5.44e-07 │                  1.16e-07 │                 9.61e-08 │
│   k_rna │       3.24e-07 │           6.71e-07 │                  6.60e-07 │                 6.60e-07 │
│k_pro_sat│       4.50e-13 │           4.82e-14 │                  3.21e-11 │                 2.73e-12 │
│k_dna_sat│       8.29e-16 │           2.07e-16 │                  1.45e-12 │                 6.01e-14 │
│k_rna_sat│       1.00e-12 │           2.40e-14 │                  1.76e-10 │                 1.74e-11 │
└─────────┴────────────────┴────────────────────┴───────────────────────────┴──────────────────────────┘
=#
"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::MacroMolecularMode)
    params=Dict(
        "κhP"      => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"      => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"   => [1],       # Number of phyto cells each super individual represents
        "Cquota"   => [1.8e-13], # C quota of phyto cells (mmolC/cell)
        "C_DNA"    => [1.8e-13], # DNA C quota of phyto cells (mmolC/cell)
        "var"      => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "RNA2DNA"  => [1.73],    # Initial RNA:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "PRO2DNA"  => [34.83],   # Initial protein:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "CH2DNA"   => [23.3],    # Initial carbohydrate+lipid:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "Chl2DNA"  => [3.45],    # Initial Chla:DNA ratio in phytoplankton (mmolC/mmolC) from Micromonas sp.
        "α"        => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"        => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"     => [27.0],    # Optimal temperature for growth (C)
        "Tmax"     => [30.0],    # Maximal temperature for growth (C)
        "Ea"       => [5.3e4],   # Free energy
        "PCmax"    => [6.2e-5],  # Maximum primary production rate (per second)
        "VDOCmax"  => [0.0],     # Maximum DOC uptake rate (mmol C/mmol C/second)
        "VNH4max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VNO3max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VPO4max"  => [1.2e-6],  # Maximum P uptake rate (mmol P/mmol C/second)
        "KsatNH4"  => [0.005],   # Half-saturation coeff (mmol N/m³)
        "KsatNO3"  => [0.010],   # Half-saturation coeff (mmol N/m³)
        "KsatPO4"  => [0.003],   # Half-saturation coeff (mmol P/m³)
        "KsatDOC"  => [0.0],     # Half-saturation coeff (mmol C/m³)
        "NSTmax"   => [0.12],    # Maximum N reserve in total N (mmol N/mmol N)
        "PSTmax"   => [0.70],    # Maximum P reserve in total P (mmol P/mmol P)
        "CHmax"    => [0.4],     # Maximum Carbohydrate in cell (mmol C/mmol C)
        "respir"   => [1.2e-6],  # Respiration rate(per second)
        "k_pro"    => [6.0e-5],  # Protein synthesis rate (mmol C/mmol C/second)
        "k_sat_pro"=> [4.5e-13], # Hafl saturation constent for protein synthesis (mmol C/cell)
        "k_dna"    => [1.6e-7],  # DNA synthesis rate (mmol C/mmol C/second)
        "k_sat_dna"=> [1.0e-15], # Hafl saturation constent for DNA synthesis (mmol C/cell)
        "k_rna"    => [3.0e-7],  # RNA synthesis rate (mmol C/mmol C/second)
        "k_sat_rna"=> [1.0e-12], # Hafl saturation constent for RNA synthesis (mmol C/cell)
        "Chl2N"    => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC_PRO" => [1/4.5],   # N:C ratio in protein (from Inomura et al 2020.)
        "R_NC_DNA" => [1/2.9],   # N:C ratio in DNA (from Inomura et al 2020.)
        "R_PC_DNA" => [1/11.1],  # P:C ratio in DNA
        "R_NC_RNA" => [1/2.8],   # N:C ratio in RNA (from Inomura et al 2020.)
        "R_PC_RNA" => [1/10.7],  # P:C ratio in RNA
        "dvid_P"   => [1.0e-5],  # Division probability per second
        "dvid_reg" => [2.0],     # Regulation of cell division
        "grz_P"    => [0.0],     # Grazing probability per second
        "mort_P"   => [5e-5],    # Probability of cell natural death per second
        "mort_reg" => [0.5],     # Regulation of cell natural death
        "grazFracC"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP"=> [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP"=> [0.5],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end


"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::IronEnergyMode)
    params=Dict(
        "κhP"       => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"       => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"    => [1],       # Number of phyto cells each super individual represents (cells)
        "Cquota"    => [1.8e-11], # C quota of phyto cells at size = 1.0 (mmolC/cell)
        "Rad"       => [0.12],    # Radius (μm) of Prochlorococcus
        "mean"      => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"       => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2Cint"  => [0.10],    # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"         => [4.5e-2],  # Irradiance absorption coeff (mmolC m² second/mgChl /μmol photon)
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea"        => [5.3e4],   # Free energy
        "is_nr"     => [1.0],     # 1 for non-diazotroph, 0 for diazotroph
        "is_croc"   => [0.0],     # 1 for Crocosphaera-like N fixation pattern
        "is_tric"   => [0.0],     # 1 for Trichodesmium-like N fixation pattern
        "PCmax"     => [8.0e-5],  # Maximum light harvesting rate (mmolATP/mmolC/second)
        "VNH4max"   => [6.9e-6],  # Maximum N uptake rate (mmolN/mmolC/second)
        "VNO3max"   => [6.9e-6],  # Maximum N uptake rate (mmolN/mmolC/second)
        "VPO4max"   => [1.2e-6],  # Maximum P uptake rate (mmolP/mmolC/second)
        "k_O2"      => [3.0e-12], # O₂ permeability coefficient (m²/s)
        "k_cf"      => [1.0e-5],  # Carbon fixation rate (per second)
        "k_rs"      => [1.5e-6],  # Maximum respiration rate (per second)
        "k_nr"      => [1.5e-5],  # Nitrate reduction rate (per second)
        "k_nf"      => [2.8e-6],  # N fixation rate (mmolN/mmolC/second)
        "k_mtb"     => [3.5e-5],  # Metabolic rate (per second)
        "e_cf"      => [3.0],     # Energy consumption ratio of carbon fixation (mmolATP/mmolC)
        "e_rs"      => [5.0],     # Energy production ratio of respiration (mmolATP/mmolC)
        "e_nf"      => [8.0],     # Energy consumption ratio of N fixation (mmolATP/mmolN)
        "e_nr"      => [0.0],     # Energy consumption ratio of NO3 reduction (mmolATP/mmolN)
        "re_ps"     => [0.77],    # NADPH production ratio of light harvesting (mmolNADPH/mmolATP)
        "re_cf"     => [2.0],     # NADPH consumption ratio of carbon fixation (mmolNADPH/mmolC)
        "re_rs"     => [0.33],    # NADPH production ratio of respiration (mmolNADPH/mmolC)
        "re_nf"     => [2.0],     # NADPH consumption ratio of N fixation (mmolNADPH/mmolN)
        "re_nr"     => [4.0],     # NADPH consumption ratio of NO3 reduction (mmolNADPH/mmolN)
        "o_ps"      => [0.38],    # O2 production ratio of light harvesting (mmolO2/mmolATP)
        "k_Fe_ST2PS"=> [2.4e-5],  # Allocation rate of Fe from storage to PS (per second)
        "k_Fe_PS2ST"=> [1.2e-6],  # Allocation rate of Fe from PS to storage (per second)
        "k_Fe_ST2NR"=> [1.2e-5],  # Allocation rate of Fe from storage to NR (per second)
        "k_Fe_NR2ST"=> [1.2e-5],  # Allocation rate of Fe from NR to storage (per second)
        "k_Fe_ST2NF"=> [1.2e-5],  # Allocation rate of Fe from storage to NF (per second)
        "k_Fe_NF2ST"=> [1.2e-5],  # Allocation rate of Fe from NF to storage (per second)
        "KfePS"     => [3.0e-6],  # Haff-saturation coeff of iron quota for photosynthesis (mmolFe/mmolC)
        "KfeNR"     => [2.0e-6],  # Haff-saturation coeff of iron quota for NO3 reduction (mmolFe/mmolC)
        "KfeNF"     => [5.0e-6],  # Haff-saturation coeff of iron quota for N fixation (mmolFe/mmolC)
        "KsatNH4"   => [0.005],   # Half-saturation coeff (mmolN/m³)
        "KsatNO3"   => [0.010],   # Half-saturation coeff (mmolN/m³)
        "KsatPO4"   => [0.003],   # Half-saturation coeff (mmolP/m³)
        "KSAFe"     => [2.77e-7], # Surface-area specific iron uptake rate (m/cell/second)
        "qNO3max"   => [0.25],    # Maximum NO3 quota in cell (mmolN/mmolC)
        "qNH4max"   => [0.25],    # Maximum NH4 quota in cell (mmolN/mmolC)
        "qPmax"     => [0.02],    # Maximum P quota in cell (mmolP/mmolC)
        "qFemax"    => [2.0e-5],  # Maximum Fe quota in cell (mmolFe/mmolC)
        "CHmax"     => [0.4],     # Maximum C quota in cell (mmolC/mmolC)
        "qO2diff"   => [2.0e2],   # Intracellular O2 concentration when O2 diffusion reaches maximum (mmolO₂/m³)
        "qO2nf"     => [1.0e2],   # Intracellular O2 concentration when N fixation reaches 0.0 (mmolO₂/m³)
        "Chl2N"     => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC"      => [16/106],  # N:C ratio in cell biomass
        "R_PC"      => [1/106],   # N:C ratio in cell biomass
        "NF_clock"  => [21600.0], # the circadian clock for N fixation
        "grz_P"     => [0.0],     # Grazing probability per second
        "dvid_P"    => [1e-4],    # Probability of cell division per second.
        "dvid_type" => [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_reg"  => [2.5],     # Regulations of cell division (cell size)
        "dvid_reg2" => [12.0],    # Regulations of cell division (clock time)
        "mort_P"    => [5e-5],    # Probability of cell natural death per second
        "mort_reg"  => [0.5],     # Regulation of cell natural death
        "grazFracC" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracFe"=> [0.1],     # Fraction goes into dissolved organic pool
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracFe"=> [0.1],     # Fraction goes into dissolved organic pool
        "ICPE_ptc"  => [5e-5],    # Photoelectrochemical conversion efficiency of Photoelectric Particles
        "SA_e"      => [4.5e-14], # Projected area (m²) Prochlorococcus
        "eATP"      => [7.5e-5],  # Efficiency of electron-to-ATP conversion (mmolATP/µmol electron)
        "sz_min"    => [1.0e-3],  # Minimal size of a abiotic particle (mmolFe/particle)
        "ptc_de"    => [4560.0],  # Density of iron mineral (kg/m³)
        "Fe_frac"   => [0.65],    # Mass fraction of Fe in iron mineral
        "M_Fe"      => [55.85],   # Molar mass of Fe (g/mol)
        "max_ptc"   => [25],      # maximum number of abiotic particles that can interact with one phytoplankton cell)
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end


"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::QuotaMode)
    params=Dict(
        "κhP"      => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"      => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"   => [1],       # Number of phyto cells each super individual represents
        "Cquota"   => [1.8e-11], # C quota of phyto cells at size = 1.0
        "mean"     => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"      => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2Cint" => [0.10],    # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"        => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"        => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"     => [27.0],    # Optimal temperature for growth (C)
        "Tmax"     => [30.0],    # Maximal temperature for growth (C)
        "Ea"       => [5.3e4],   # Free energy
        "PCmax"    => [4.2e-5],  # Maximum primary production rate (per second)
        "VDOCmax"  => [0.0],     # Maximum DOC uptake rate (mmol C/mmol C/second)
        "VNH4max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VNO3max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VPO4max"  => [1.2e-6],  # Maximum P uptake rate (mmol P/mmol C/second)
        "KsatNH4"  => [0.005],   # Half-saturation coeff (mmol N/m³)
        "KsatNO3"  => [0.010],   # Half-saturation coeff (mmol N/m³)
        "KsatPO4"  => [0.003],   # Half-saturation coeff (mmol P/m³)
        "KsatDOC"  => [0.0],     # Half-saturation coeff (mmol C/m³)
        "Nqmax"    => [0.25],    # Maximum N quota in cell (mmol N/mmol C)
        "Pqmax"    => [0.02],    # Maximum P quota in cell (mmol P/mmol C)
        "Cqmax"    => [0.4],     # Maximum C quota in cell (mmol C/mmol C)
        "k_mtb"    => [3.5e-5],  # Metabolic rate (per second)
        "respir"   => [1.2e-6],  # Respiration rate(per second)
        "Chl2N"    => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC"     => [16/106],  # N:C ratio in cell biomass
        "R_PC"     => [1/106],   # N:C ratio in cell biomass
        "grz_P"    => [0.0],     # Grazing probability per second
        "dvid_P"   => [1e-4],    # Probability of cell division per second.
        "dvid_type"=> [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_reg" => [2.5],     # Regulations of cell division (cell size)
        "dvid_reg2"=> [12.0],    # Regulations of cell division (clock time)
        "mort_P"   => [5e-5],    # Probability of cell natural death per second
        "mort_reg" => [0.5],     # Regulation of cell natural death
        "grazFracC"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP"=> [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP"=> [0.5],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end


"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::CarbonMode)
    params=Dict(
        "κhP"       => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"       => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"    => [1],       # Number of phyto cells each super individual represents
        "Cquota"    => [1.8e-11], # C quota of phyto cells at size = 1.0
        "mean"      => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"       => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2C"     => [0.10],    # Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"         => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"         => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea_ref"    => [5.3e4],   # Free energy
        "Topt_ref"  => [26.0],    # reference optimal growth temperature(C)
        "PCmax"     => [4.2e-5],  # Maximum primary production rate (per second)
        "respir"    => [1.2e-6],  # Respiration rate(per second)
        "f_T2B"     => [2.7e-6],  # Thermal damage rate (1.0/K/s)
        "grz_P"     => [0.0],     # Grazing probability per second
        "dvid_P"    => [5e-5],    # Probability of cell division per second.
        "dvid_type" => [1],       # The type of cell division, 1:sizer, 2:adder
        "dvid_reg"  => [1.9],     # Regulations of cell division (sizer)
        "dvid_reg2" => [12.0],    # Regulations of cell division (sizer)
        "mort_P"    => [5e-5],    # Probability of cell natural death per second
        "mort_reg"  => [5.0e-2],  # Regulation of cell natural death
        "grazFracC" => [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
        "thermal"   => [1.0],     # thermal damage, 1 for on, 0 for off
        "is_bact"   => [0.0],     # is this specis bacteria or not, 1 for bacteria, 0 for phytoplankton
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end


"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::ProteinMode)
    params=Dict(
        "κhP"        => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"        => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"     => [1],       # Number of phyto cells each super individual represents
        "Cquota"     => [1.8e-11], # Structural C quota of phyto cells at size = 1.0 (mmolC/cell)
        "C_DNA"      => [1.8e-13], # DNA C quota of phyto cells (mmolC/cell)
        "var"        => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "RNA2DNA"    => [1.73],    # Initial RNA:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "PRO_RB2DNA" => [2.84],    # Initial ribosomal protein:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_MC2DNA" => [1.71],    # Initial carbon metabolic protein:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_MN2DNA" => [2.9e-2],  # Initial nitrogen metabolism protein:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_TN2DNA" => [8.9e-2],  # Initial nitrate transporter:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_TP2DNA" => [5.6e-2],  # Initial phosphate transporter:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_TFe2DNA"=> [8.5e-4],  # Initial iron transporter:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_RS2DNA" => [2.91],    # Initial carbon metabolism protein:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_PS2DNA" => [7.03],    # Initial photosynthesis protein:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "PRO_OT2DNA" => [293.94],  # Initial other protein:DNA ratio in phytoplankton (mmol C/mmolC) from Synechococcus sp.
        "CH2DNA"     => [23.3],    # Initial carbohydrate+lipid:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "Chl2DNA"    => [3.45],    # Initial Chla:DNA ratio in phytoplankton (mmolC/mmolC) from Micromonas sp.
        "α"          => [8e-6],    # Irradiance absorption coeff (mmolATP m² /mgChl /μmol photon)
        "Topt"       => [27.0],    # Optimal temperature for growth (C)
        "Tmax"       => [30.0],    # Maximal temperature for growth (C)
        "Ea"         => [5.3e4],   # Free energy
        "is_nr"      => [1.0],     # 1 for non-diazotroph, 0 for diazotroph
        "is_croc"    => [0.0],     # 1 for Crocosphaera-like N fixation pattern
        "is_tric"    => [0.0],     # 1 for Trichodesmium-like N fixation pattern
        "PCmax"      => [7.8e-8],  # Maximum light harvesting rate (mmolATP/mmolC/second)
        "VDOCmax"    => [0.0],     # Maximum DOC uptake rate (mmol C/mmol C/second)
        "KcatCF"     => [7.3e-4],  # Maximum carbon fixation rate per carbon fixation protein (per second) from Synechococcus sp.
        "KcatNF"     => [1.0],     # Maximum N fixation rate per N fixation protein (mmolN/mmolC/second)
        "KcatNR"     => [1.3e-2],  # Maximum nitrate reduction rate per nitrate reduction protein (mmolN/mmolC/second) from Synechococcus sp.
        "KcatTNH4"   => [1.9e-4],  # Maximum ammonium uptake rate per NH4 transporter protein (mmolN/mmolC/second) from bacteria
        "KcatTNO3"   => [2.4e-4],  # Maximum nitrate uptake rate per NO3 transporter protein (mmolN/mmolC/second) from bacteria
        "KcatTPO4"   => [2.7e-4],  # Maximum phosphate uptake rate per PO4 transporter protein (mmolP/mmolC/second) from bacteria
        "KcatTFe"    => [8.5e-5],  # Maximum iron uptake rate per Fe transporter protein (mmolFe/mmolC/second) from Synechococcus sp.
        "KcatRS"     => [3.0e-5],  # Maximum respiration rate per respiration protein (per second) from Synechococcus sp.
        "KcatRB"     => [3.0e-3],  # Maximum ribosome synthesis rate per ribosome protein (per second) from E. coli
        "β_MC"       => [0.12],    # Fraction of C metabolic enzyme complex per C metabolic protein from Synechococcus sp.
        "β_MN"       => [2e-3],    # Fraction of ribosomes synthesizing carbon metabolic enzyme complex  from Synechococcus sp.
        "β_RB"       => [0.19],    # Fraction of ribosomes synthesizing ribosomal protein from Synechococcus sp.
        "β_TN"       => [6e-3],    # Fraction of ribosomes synthesizing ammonium transporter protein from Synechococcus sp.
        "β_TPO4"     => [4e-3],    # Fraction of ribosomes synthesizing phosphate transporter protein from Synechococcus sp.
        "β_TFe"      => [6e-5],    # Fraction of ribosomes synthesizing iron transporter protein from Synechococcus sp.
        "β_RS"       => [0.2],     # Fraction of ribosomes synthesizing respiration protein from Synechococcus sp.
        "β_PS"       => [0.48],    # Fraction of ribosomes synthesizing photosynthesis protein from Synechococcus sp.
        "β_OT"       => [0.05],    # Fraction of ribosomes synthesizing other proteins from Synechococcus sp.
        "KsatNH4"    => [0.005],   # Half-saturation coeff (mmol N/m³)
        "KsatNO3"    => [0.010],   # Half-saturation coeff (mmol N/m³)
        "KsatPO4"    => [0.003],   # Half-saturation coeff (mmol P/m³)
        "KsatDOC"    => [0.0],     # Half-saturation coeff (mmol C/m³)
        "KsatNR"     => [0.005],   # Half-saturation coeff (mmol N / mmolC)
        "KsatFe"     => [1.0e-6],  # Half-saturation coeff of iron quota for iron transport (mmolFe/m³)
        "KFe_em"     => [1.0e-5],  # Half-saturation coeff of iron quota for photosynthesis, nitrogen fixation, nitrate reduction and respiration (mmolFe/mmolC)
        "TN_max"     => [200],     # Maximum nitrate concentration for protein synthesis in cell (mmol N/m3)
        "TP_max"     => [20],      # Maximum phosphate concentration for protein synthesis in cell (mmol P/m3)
        "TFe_max"    => [1.0e-2],  # Maximum iron concentration for protein synthesis in cell (mmol Fe/m3)
        "PSTmax"     => [0.05],    # Maximum P reserve in cell (mmol P/mmol C)
        "CHmax"      => [0.4],     # Maximum Carbohydrate in cell (mmol C/mmol C)
        "qNH4max"    => [0.25],    # Maximum NH4 quota in cell (mmolN/mmolC)
        "qNO3max"    => [0.25],    # Maximum NO3 quota in cell (mmolN/mmolC)
        "qFemax"     => [2.0e-5],  # Maximum Fe quota in cell (mmolFe/mmolC)
        "PARmax"     => [600.0],   # Maximum PAR for photosynthesis protein synthesis (μmol photon/m²/second)
        "QC2N_dna"   => [10],       # Maximum cellular carbon to nitrogen ratio (mmolC/mmolN)
        "QC2N_chl"   => [10],      # Maximum cellular carbon to nitrogen ratio (mmolC/mmolN)
        "QN2C_DP"    => [0.25],    # Maximum cellular nitrogen to carbon ratio for DP (mmolN/mmolC)
        "k_degRB"    => [2.8e-5],  # Ribosome protein degradation rate (per second)
        "k_degPS"    => [1.9e-5],  # photosynthesis protein degradation rate (per second)
        "k_degMC"    => [1.1e-5],  # Carbon fixation protein degradation rate (per second)
        "k_degChl"   => [3.3e-5],  # Chlorophyll degradation rate (per second)
        "k_degRNA"   => [2.8e-5],  # RNA degradation rate (per second)
        "lag_shape"  => [10.8],    # Shape parameter for the lag distribution (Gamma)
        "lag_scale"  => [1.89 * 86400.0],     # Scale parameter for the lag distribution (Gamma)
        "PRO_RBmin"  => [5.9e-13], # Minimum ribosome protein quota (mmol C/individual) from Synechococcus sp.
        "PRO_PSmin"  => [1.5e-12], # Minimum photosynthesis protein quota (mmol C/individual) from Synechococcus sp.
        "PRO_MCmin"  => [3.6e-13], # Minimum carbon fixation protein quota (mmol C/individual) from Synechococcus sp.
        "RNAmin"     => [7.9e-13], # Minimum RNA quota (mmol C/individual) from Synechococcus sp.
        "Chlmin"     => [9e-12],   # Minimum Chl quota (mg Chl/individual)
        "e_TNO3"     => [1.0],     # Energy consumption rate of nitrate transporter (mmolATP/mmolN) from Synechococcus sp.
        "e_TPO4"     => [1.0],     # Energy consumption rate of phosphate transporter (mmolATP/mmolP) from Synechococcus sp.
        "e_TFe"      => [1.0],     # Energy consumption rate of iron transporter (mmolATP/mmolFe) from Synechococcus sp.
        "e_CF"       => [8.0],     # Energy consumption rate of carbon fixation (mmolATP/mmolC)
        "e_RS"       => [5.0],     # Energy production rate of respiration (mmolATP/mmolC)
        "e_NF"       => [8.0],     # Energy consumption rate of N fixation (mmolATP/mmolN)
        "e_NR"       => [10.0],    # Energy consumption rate of NO3 reduction (mmolATP/mmolN)
        "e_SP"       => [1.55],    # Energy consumption rate of protein synthesis (mmolATP/mmol C) from Synechococcus sp.
        "e_DNA"      => [1.17],    # Energy consumption rate of DNA synthesis (mmolATP/mmol C) from Synechococcus sp.
        "e_RNA"      => [0.98],    # Energy consumption rate of RNA synthesis (mmolATP/mmol C) from Synechococcus sp.
        "e_min"      => [7.72e-17],# Minimum maintenance energy requirement (mmolATP/individual/second) from E.coli
        "k_DNA"      => [1.6e-7],  # DNA synthesis rate (mmol C/second)
        "k_sat_DNA"  => [1.0e-15], # Half-saturation constant for DNA synthesis (mmol C/cell)
        "k_sat_RNA"  => [1.0e-12], # Half-saturation constant for RNA synthesis (mmol C/cell)
        "k_sat_PRO"  => [4.5e-13], # Half-saturation constant for protein synthesis (mmol C/cell)
        "Chl2C"      => [0.54],    # Maximum Chla:C ratio in phytoplankton (mgChl/mmolC)
        "R_NC_PRO"   => [1/4.5],   # N:C ratio in protein (from Inomura et al 2020.)
        "R_NC_DNA"   => [1/2.9],   # N:C ratio in DNA (from Inomura et al 2020.)
        "R_PC_DNA"   => [1/11.1],  # P:C ratio in DNA
        "R_NC_RNA"   => [1/2.8],   # N:C ratio in RNA (from Inomura et al 2020.)
        "R_PC_RNA"   => [1/10.7],  # P:C ratio in RNA
        "R_C_RNAPRB" => [1.53],    # C:C ratio in RNA/ ribosome protein
        "dvid_P"     => [1.0e-5],  # Division probability per second
        "dvid_reg"   => [2],       # Regulations of cell division (cell size)
        "grz_P"      => [0.0],     # Grazing probability per second
        "mort_P"     => [5e-5],    # Probability of cell natural death per second
        "mort_reg"   => [0.5],     # Regulation of cell natural death
        "grazFracC"  => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN"  => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP"  => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracFe" => [0.1],     # Fraction goes into dissolved organic pool
        "mortFracC"  => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN"  => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP"  => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracFe" => [0.1],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end


"""
    colony_params_default(Ncl::Int64, Nsp::AbstractArray, mode::AbstractMode)
Generate default colony parameter values based on `AbstractMode`, colony number `Ncl` and species number `Nsp`.
"""
function colony_params_default(Ncl::Int64, Nsp::AbstractArray, mode::IronEnergyMode)
    params=Dict(
        "κhP"       => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"       => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"    => [1],       # Number of phyto cells each super individual represents (cells)
        "Cquota"    => [1.8e-11], # C quota of phyto cells at size = 1.0 (mmolC/cell)
        "Rad"       => [0.12],    # Radius (μm) of Prochlorococcus
        "mean"      => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"       => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2Cint"  => [0.10],    # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"         => [4.5e-2],  # Irradiance absorption coeff (mmolC m² second/mgChl /μmol photon)
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea"        => [5.3e4],   # Free energy
        "is_nr"     => [1.0],     # 1 for non-diazotroph, 0 for diazotroph
        "is_croc"   => [0.0],     # 1 for Crocosphaera-like N fixation pattern
        "is_tric"   => [0.0],     # 1 for Trichodesmium-like N fixation pattern
        "PCmax"     => [8.0e-8],  # Maximum light harvesting rate (mmolATP/mmolC/second)
        "VNH4max"   => [6.9e-6],  # Maximum N uptake rate (mmolN/mmolC/second)
        "VNO3max"   => [6.9e-6],  # Maximum N uptake rate (mmolN/mmolC/second)
        "VPO4max"   => [1.2e-6],  # Maximum P uptake rate (mmolP/mmolC/second)
        "k_O2"      => [3.0e-12], # O₂ permeability coefficient (m²/s)
        "k_cf"      => [1.0e-5],  # Carbon fixation rate (per second)
        "k_rs"      => [1.5e-6],  # Maximum respiration rate (per second)
        "k_nr"      => [2.8e-6],  # Nitrate reduction rate (per second)
        "k_nf"      => [2.8e-6],  # N fixation rate (mmolN/mmolC/second)
        "k_mtb"     => [3.5e-5],  # Metabolic rate (per second)
        "e_cf"      => [3.0],     # Energy consumption ratio of carbon fixation (mmolATP/mmolC)
        "e_rs"      => [5.0],     # Energy production ratio of respiration (mmolATP/mmolC)
        "e_nf"      => [8.0],     # Energy consumption ratio of N fixation (mmolATP/mmolN)
        "e_nr"      => [0.0],     # Energy consumption ratio of NO3 reduction (mmolATP/mmolN)
        "re_ps"     => [0.76],    # NADPH production ratio of light harvesting (mmolNADPH/mmolATP)
        "re_cf"     => [2.0],     # NADPH consumption ratio of carbon fixation (mmolNADPH/mmolC)
        "re_rs"     => [0.33],    # NADPH production ratio of respiration (mmolNADPH/mmolC)
        "re_nf"     => [2.0],     # NADPH consumption ratio of N fixation (mmolNADPH/mmolN)
        "re_nr"     => [4.0],     # NADPH consumption ratio of NO3 reduction (mmolNADPH/mmolN)
        "o_ps"      => [0.38],    # O2 production ratio of light harvesting (mmolO2/mmolATP)
        "k_Fe_ST2PS"=> [2.4e-5],  # Allocation rate of Fe from storage to PS (per second)
        "k_Fe_PS2ST"=> [1.2e-6],  # Allocation rate of Fe from PS to storage (per second)
        "k_Fe_ST2NR"=> [1.2e-5],  # Allocation rate of Fe from storage to NR (per second)
        "k_Fe_NR2ST"=> [1.2e-5],  # Allocation rate of Fe from NR to storage (per second)
        "k_Fe_ST2NF"=> [1.2e-5],  # Allocation rate of Fe from storage to NF (per second)
        "k_Fe_NF2ST"=> [1.2e-5],  # Allocation rate of Fe from NF to storage (per second)
        "KfePS"     => [3.0e-6],  # Haff-saturation coeff of iron quota for photosynthesis (mmolFe/mmolC)
        "KfeNR"     => [2.0e-6],  # Haff-saturation coeff of iron quota for NO3 reduction (mmolFe/mmolC)
        "KfeNF"     => [5.0e-6],  # Haff-saturation coeff of iron quota for N fixation (mmolFe/mmolC)
        "KsatNH4"   => [0.005],   # Half-saturation coeff (mmolN/m³)
        "KsatNO3"   => [0.010],   # Half-saturation coeff (mmolN/m³)
        "KsatPO4"   => [0.003],   # Half-saturation coeff (mmolP/m³)
        "KSAFe"     => [2.77e-7], # Surface-area specific iron uptake rate (m/cell/second)
        "qNO3max"   => [0.25],    # Maximum NO3 quota in cell (mmolN/mmolC)
        "qNH4max"   => [0.25],    # Maximum NH4 quota in cell (mmolN/mmolC)
        "qPmax"     => [0.02],    # Maximum P quota in cell (mmolP/mmolC)
        "qFemax"    => [2.0e-5],  # Maximum Fe quota in cell (mmolFe/mmolC)
        "CHmax"     => [0.4],     # Maximum C quota in cell (mmolC/mmolC)
        "qO2diff"   => [2.0e2],   # Intracellular O2 concentration when O2 diffusion reaches maximum (mmolO₂/m³)
        "qO2nf"     => [1.0e2],   # Intracellular O2 concentration when N fixation reaches 0.0 (mmolO₂/m³)
        "Chl2N"     => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC"      => [16/106],  # N:C ratio in cell biomass
        "R_PC"      => [1/106],   # N:C ratio in cell biomass
        "NF_clock"  => [21600.0], # the circadian clock for N fixation
        "grz_P"     => [0.0],     # Grazing probability per second
        "dvid_P"    => [1e-4],    # Probability of cell division per second.
        "dvid_type" => [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_reg"  => [2.5],     # Regulations of cell division (cell size)
        "dvid_reg2" => [12.0],    # Regulations of cell division (clock time)
        "mort_P"    => [5e-5],    # Probability of cell natural death per second
        "mort_reg"  => [0.5],     # Regulation of cell natural death
        "grazFracC" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracFe"=> [0.1],     # Fraction goes into dissolved organic pool
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracFe"=> [0.1],     # Fraction goes into dissolved organic pool
        "k_TqP"     => [0.05],    # Phosphorus exchange rate between cells within a colony (per second)
        "k_TCH"     => [0.05],    # CH exchange rate between cells within a colony (per second)
        "k_TqO2"    => [0.05],    # O₂ exchange rate between cells within a colony (per second)
        "k_TqFe"    => [0.05],    # iron exchange rate between cells within a colony (per second)
        "k_TqNH4"   => [0.05],    # NH4 exchange rate between cells within a colony (per second)
        "ϵ_TqP"     => [0.25],    # Phosphorus exchange cost between cells within a colony (per second)
        "ϵ_TCH"     => [0.25],    # CH exchange cost between cells within a colony (per second)
        "ϵ_TqO2"    => [0.25],    # O₂ exchange cost between cells within a colony (per second)
        "ϵ_TqFe"    => [0.25],    # iron exchange cost between cells within a colony (per second)
        "ϵ_TqNH4"   => [0.25],    # NH4 exchange cost between cells within a colony (per second)
    )

    param_cl = []
    for i in 1:Ncl
        if Nsp[i] > 1
            params = generate_n_species_params(Nsp[i], params)
        end
        push!(param_cl, copy(params))
    end
    return param_cl

end


"""
    abiotic_params_default(N::Int64)
Generate default abiotic particle parameter values based on species number `N`.
"""
function abiotic_params_default(N::Int64)
    params=Dict(
        "κhP"       => [0.0],     # Horizontal particle diffusivity (m²/s)
        "κvP"       => [0.0],     # Vertical particle diffusivity (m²/s)
        "Nsuper"    => [1],       # Number of abiotic particles each super individual represents
        "Rd"        => [1.0e-4],  # Distance between abiotic particle and phytoplankton cell (m)
        "release_P" => [1.0e-6],  # Probability of particle release per second
        "sz_min"    => [1.0e-3],  # Minimal size of a abiotic particle (mmolFe/particle)
        "Ktr"       => [0.0e-5],  # Rate of particle forming from tracer
    )
    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    update_phyt_params(tmp::Dict, FT::DataType; N::Int64, mode::AbstractMode)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
- `mode` is the mode of phytoplankton physiology resolved in the model
"""
function update_phyt_params(tmp::Dict, FT::DataType; N::Int = 1, mode::AbstractMode = QuotaMode())
    parameters = phyt_params_default(N,mode)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: phyt parameter not found $key"))
        else
            parameters[key] = FT.(tmp[key])
        end
    end
    return parameters
end

"""
    update_colony_params(tmp::Dict, FT::DataType; N::Int64, mode::AbstractMode)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
- `mode` is the mode of phytoplankton physiology resolved in the model
"""
function update_colony_params(tmps::AbstractArray, FT::DataType;
                              Ncl::Int = 1, Nsp::AbstractArray = [1],
                              mode::AbstractMode = IronEnergyMode())
    parameters = colony_params_default(Ncl, Nsp, mode)
    for i in eachindex(tmps)
        tmp = tmps[i]
        tmp_keys = collect(keys(tmp))
        pkeys = collect(keys(parameters[i]))
        for key in tmp_keys
            if length(findall(x->x==key, pkeys))==0
                throw(ArgumentError("PARAM: colony parameter not found $key"))
            else
                parameters[i][key] = FT.(tmp[key])
            end
        end
    end
    return parameters
end


"""
    update_abiotic_params(tmp::Dict, FT::DataType; N::Int64)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
"""
function update_abiotic_params(tmp::Dict, FT::DataType; N::Int = 1)
    parameters = abiotic_params_default(N)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: abiotic parameter not found $key"))
        else
            parameters[key] = FT.(tmp[key])
        end
    end
    return parameters
end
