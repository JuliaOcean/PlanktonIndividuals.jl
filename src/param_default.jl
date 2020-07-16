# Parameters
param_default=Dict(
    "P_Nind"   => 1000,               # Number of phyto individuals of each species
    "P_Nsp"    => 1,                  # Number of phyto species
    "P_Nsuper" => 1e10,               # Number of phyto cells each super individual represents
    "P_Cquota" => [1.8e-11, 1.8e-10], # C quota of phyto cells at size = 1.0
    "useDOC"   => [0,0],              # A switch that controls whether the individual can uptake DOC
    "isDiaz"   => [0,0],              # A switch that controls whether the individual can fix N
    "P_mean"   => 1.5,                # Mean of the normal distribution of initial phyto individuals
    "P_var"    => 0.3,                # Variance of the normal distribution of initial phyto individuals
    "Z_Nind"   => 0,                  # Number of zoo individuals of each species
    "Z_Nsp"    => 1,                  # Number of zoo species
    "Z_Nsuper" => 1e0,                # Number of zoo cells each super individual represents
    "Z_Cquota" => 0.0,                # C quota of zoo cells at size = 1.0
    "Z_mean"   => 0.0,                # Mean of the normal distribution of initial zoo individuals
    "Z_var"    => 0.0,                # Variance of the normal distribution of initial zoo individuals
    "diag_inds"=> zeros(Int,14),      # Diagnostic indices, refer to diagnostics.jl
    "diag_freq"=> 3600,               # Frequency of diagnostics (second)
    "PCmax"    => [1.8, 1.8],         # Maximum primary production rate (per day)
    "PC_b"     => [0.6, 0.6],         # Shape parameter for size
    "Chl2N"    => 3.0,                # Maximum Chla:N ratio in phytoplankton
    "Chl2Cint" => [0.16,0.16],        # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
    "R_NC"     => 16/106,             # N:C ratio in cell biomass
    "R_PC"     => 1/106,              # N:C ratio in cell biomass
    "α"        => [2.0e-2,2.0e-2],    # Irradiance absorption coeff (m²/mgChl)
    "Φ"        => [4.0e-5,4.0e-5],    # Maximum quantum yield (mmolC/μmol photon)
    "katten_w" => 0.046,              # PAR attenuation (/m)
    "katten_c" => 0.04,               # PAR attenuation (m²/mgChl)
    "inhibcoef"=> [0.0, 0.0],         # Light inhibition
    "Grz_P"    => 4000,               # Grazing probability scalar
    "P_dvid"   => [0.2,0.2],          # Probability scaler of cell division.
    "dvid_type"=> [1,1],              # The type of cell division, 1:sizer, 2:adder.
    "dvid_stp" => [6.0,6.0],          # Steepness of sigmoidal function
    "dvid_reg" => [1.9,1.9],          # Regulations of cell division(sizer, adder, timer)
    "P_death"  => [0.3,0.3],          # Probability scaler of cell natural death
    "death_reg"=> [0.5,0.5],          # Regulation of cell natural death
    "Tempref"  => 293.15,             # Reference temperature in K
    "TempAe"   => -4000.0,            # Arrenhius equation
    "TempCoeff"=> 0.8,                # Arrenhius equation
    "VNH4max"  => [0.6, 0.6],         # Maximum N uptake rate (mmol N/mmol C/day)
    "VNO3max"  => [0.6, 0.6],         # Maximum N uptake rate (mmol N/mmol C/day)
    "VPmax"    => [0.1, 0.1],         # Maximum P uptake rate (mmol P/mmol C/day)
    "VDOCmax"  => [0.0, 1.0],         # Maximum DOC uptake rate (mmol C/mmol C/day)
    "C2Nfix"   => [1e-2,1e-2],        # Energy cost by N fixation (mmol C/mmol N)
    "Nqmax"    => [0.17,0.17],        # Maximum N quota in cell (mmol N/mmol C)
    "Nqmin"    => [0.05,0.05],        # Minimum N quota in cell (mmol N/mmol C)
    "Pqmax"    => [0.02,0.02],        # Maximum P quota in cell (mmol P/mmol C)
    "Pqmin"    => [0.004,0.004],      # Minimum P quota in cell (mmol P/mmol C)
    "Cqmax"    => [0.4,0.4],          # Maximum C quota in cell (mmol C/mmol C)
    "Cqmin"    => [0.1,0.1],          # Minimum C quota in cell (mmol C/mmol C)
    "KsatNH4"  => [0.005, 0.005],     # Half-saturation coeff (mmol N/m³)
    "KsatNO3"  => [0.010, 0.010],     # Half-saturation coeff (mmol N/m³)
    "KsatP"    => [0.003, 0.003],     # Half-saturation coeff (mmol P/m³)
    "KsatDOC"  => [0.0, 0.3],         # Half-saturation coeff (mmol C/m³)
    "VN_b"     => [0.6, 0.6],         # Shape parameter for size
    "VP_b"     => [0.6, 0.6],         # Shape parameter for size
    "VDOC_b"   => [0.6, 0.6],         # Shape parameter for size
    "k_mtb"    => [1.0,1.0],          # Metabolic rate (per day)
    "b_k_mtb"  => [0.5,0.5],          # Metabolic rate
    "respir_a" => [0.3,0.3],          # Respiration rate(per day)
    "respir_b" => [0.13,0.13],        # Shape parameter for size
    "grazFracC"=> 0.7,                # Fraction goes into dissolved organic pool
    "grazFracN"=> 0.7,                # Fraction goes into dissolved organic pool
    "grazFracP"=> 0.7,                # Fraction goes into dissolved organic pool
    "mortFracC"=> 0.5,                # Fraction goes into dissolved organic pool
    "mortFracN"=> 0.5,                # Fraction goes into dissolved organic pool
    "mortFracP"=> 0.5,                # Fraction goes into dissolved organic pool
    "slpyFracC"=> 0.4,                # Fraction goes into zooplankton
    "slpyFracN"=> 0.4,                # Fraction goes into zooplankton
    "slpyFracP"=> 0.4,                # Fraction goes into zooplankton
    "k_sink"   => 0.01/86400,         # Sink rates for agents (m/s)
    "kDOC"     => 1/30/86400,         # Remineralization rate for DOC, turn over time: a month (per second)
    "Nit"      => 1/10/86400,         # Nitrification rate for NH4
    "kDON"     => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
    "kDOP"     => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
    "kPOC"     => 1/30/86400,         # Remineralization rate for POC, turn over time: a month (per second)
    "kPON"     => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
    "kPOP"     => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
    "κh"       => 0,                  # Horizontal diffusion
    "κv"       => 0,                  # Vertical diffusion
    "κhP"      => 0,                  # Horizontal diffusion for individuals
    "κvP"      => 0,                  # Vertical diffusion for individuals
    "v_zoo"    => 1.0e-5              # Velocity of zooplankton, m/s
)

# Options
#                    GridChoice, Gridoff, VelChoice, Veloff
RunOption=RunOptions(false,      Dict(),  false,     Dict())

#                  nTime, ΔT,   params,        Zoo
RunParam=RunParams(12,    600,  param_default, false)
