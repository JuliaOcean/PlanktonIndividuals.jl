# Parameters
function bgc_params_default()
    params = Dict(
        "kw"       => 0.046,              # PAR attenuation (/m)
        "kc"       => 0.04,               # PAR attenuation (m²/mgChl)
        "kDOC"     => 1/30/86400,         # Remineralization rate for DOC, turn over time: a month (per second)
        "Nit"      => 1/30/86400,         # Nitrification rate for NH4
        "kDON"     => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
        "kDOP"     => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
        "kPOC"     => 1/30/86400,         # Remineralization rate for POC, turn over time: a month (per second)
        "kPON"     => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
        "kPOP"     => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
        "κh"       => 1.0e-6,             # Horizontal diffusion
        "κv"       => 1.0e-6,             # Vertical diffusion
        "κhP"      => 0.0e-6,             # Horizontal diffusion for individuals
        "κvP"      => 0.0e-6,             # Vertical diffusion for individuals
    )
    return params
end

function phyt_params_default()
    params=Dict(
        "Nsuper"   => [1e1, 1e1],         # Number of phyto cells each super individual represents
        "Cquota"   => [1.8e-11, 1.8e-11], # C quota of phyto cells at size = 1.0
        "mean"     => [1.5, 1.5],         # Mean of the normal distribution of initial phyto individuals
        "var"      => [0.3, 0.3],         # Variance of the normal distribution of initial phyto individuals
        "Chl2Cint" => [0.16, 0.16],       # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"        => [2.0e-2,2.0e-2],    # Irradiance absorption coeff (m²/mgChl)
        "Φ"        => [4.0e-5,4.0e-5],    # Maximum quantum yield (mmolC/μmol photon)
        "Tempref"  => [293.15, 293.15],   # Reference temperature in K
        "TempAe"   => [-4000.0, -4000.0], # Arrenhius equation
        "TempCoeff"=> [0.8, 0.8],         # Arrenhius equation
        "PCmax"    => [2.1e-5, 2.1e-5],   # Maximum primary production rate (per second)
        "VDOCmax"  => [0.0, 1.2e-5],      # Maximum DOC uptake rate (mmol C/mmol C/second)
        "VNH4max"  => [6.9e-6, 6.9e-6],   # Maximum N uptake rate (mmol N/mmol C/second)
        "VNO3max"  => [6.9e-6, 6.9e-6],   # Maximum N uptake rate (mmol N/mmol C/second)
        "VPO4max"  => [1.2e-6, 1.2e-6],   # Maximum P uptake rate (mmol P/mmol C/second)
        "PC_b"     => [0.6, 0.6],         # Shape parameter for size
        "VDOC_b"   => [0.6, 0.6],         # Shape parameter for size
        "VN_b"     => [0.6, 0.6],         # Shape parameter for size
        "VP_b"     => [0.6, 0.6],         # Shape parameter for size
        "KsatNH4"  => [0.005, 0.005],     # Half-saturation coeff (mmol N/m³)
        "KsatNO3"  => [0.010, 0.010],     # Half-saturation coeff (mmol N/m³)
        "KsatPO4"  => [0.003, 0.003],     # Half-saturation coeff (mmol P/m³)
        "KsatDOC"  => [0.0, 0.3],         # Half-saturation coeff (mmol C/m³)
        "Nqmax"    => [0.12,0.12],        # Maximum N quota in cell (mmol N/mmol C)
        "Nqmin"    => [0.05,0.05],        # Minimum N quota in cell (mmol N/mmol C)
        "Pqmax"    => [0.01,0.01],        # Maximum P quota in cell (mmol P/mmol C)
        "Pqmin"    => [0.004,0.004],      # Minimum P quota in cell (mmol P/mmol C)
        "Cqmax"    => [0.2,0.2],          # Maximum C quota in cell (mmol C/mmol C)
        "Cqmin"    => [0.1,0.1],          # Minimum C quota in cell (mmol C/mmol C)
        "k_mtb"    => [1.2e-5,1.2e-5],    # Metabolic rate (per second)
        "k_mtb_b"  => [0.5,0.5],          # Metabolic rate
        "respir_a" => [3.5e-5,3.5e-5],    # Respiration rate(per second)
        "respir_b" => [0.13,0.13],        # Shape parameter for size
        "Chl2N"    => [3.0, 3.0],         # Maximum Chla:N ratio in phytoplankton
        "R_NC"     => [16/106, 16/106],   # N:C ratio in cell biomass
        "R_PC"     => [1/106, 1/106],     # N:C ratio in cell biomass
        "grz_P"    => [3.0e-4, 3.0e-4],   # Grazing probability scalar
        "dvid_P"   => [0.2,0.2],          # Probability scaler of cell division.
        "dvid_type"=> [1,5],              # The type of cell division, 1:sizer, 2:adder.
        "dvid_stp" => [6.0,6.0],          # Steepness of sigmoidal function
        "dvid_reg" => [1.9,1.9],          # Regulations of cell division (sizer)
        "dvid_stp2"=> [2.0,2.0],          # Steepness of sigmoidal function
        "dvid_reg2"=> [12.0,12.0],        # Regulations of cell division (sizer)
        "mort_P"   => [0.3,0.3],          # Probability scaler of cell natural death
        "mort_reg" => [0.5,0.5],          # Regulation of cell natural death
        "grazFracC"=> [0.7, 0.7],         # Fraction goes into dissolved organic pool
        "grazFracN"=> [0.7, 0.7],         # Fraction goes into dissolved organic pool
        "grazFracP"=> [0.7, 0.7],         # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5, 0.5],         # Fraction goes into dissolved organic pool
        "mortFracN"=> [0.5, 0.5],         # Fraction goes into dissolved organic pool
        "mortFracP"=> [0.5, 0.5],         # Fraction goes into dissolved organic pool
    )
    return params
end
