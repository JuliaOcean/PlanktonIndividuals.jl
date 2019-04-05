####################################### Parameters ########################################
PCmax = [2.7, 2.5] ./86400   # Maximum primary production rate (per second)
Chl2N = 3.0    # Maximum Chla:N ratio in phytoplankton
R_NC  = 16/106    # N:C ratio in cell biomass, should be lower than 16:106
Cmin  = 7.0e-13 # Minimum C quota (mmolC)
dvdcount = 0
α = 0.030       # Irradiance absorption coeff (m^2/gChl)
katten_w = 0.046 # PAR attenuation (/m)
katten_c = 0.04 # PAR attenuation (/cell)
Cquota = [1.8e-11, 1.8e-10] # Average C quota in cell (mmolC)
Grz_P  = 1500     # phyt.size/Grz_P is the probability to be grazed
dvid_size = 1.6 # relative cell size a cell can start divide
Dvid_P = 6.0   # should be greater than 2.0,phyt.size/Dvid_P is the probability to divide

Tempref = 293.15   # reference temperature in K
TempAe = -4000.0   # Arrenhius equation
TempCoeff = 0.8    # Arrenhius equation

VNmax = [1.5, 1.3] ./ 86400 # Maximum N uptake rate (per second)
Nqmax_a = 0.5   # Maximum N quota in cell (mmol N/mmol C)
Nqmin_a = 0.13  # Minimum N quota in cell (mmol N/mmol C)
KsatN = 5.0e-5  # Half-saturation coeff

respir_a = 1.0e-18 # Respiration ratio 
respir_b = 0.93    # Shape parameter

grazFracC = 0.7 # Fraction goes into dissolved organic pool
grazFracN = 0.7 # Fraction goes into dissolved organic pool
mortFracC = 0.5 # Fraction goes into dissolved organic pool
mortFracN = 0.5 # Fraction goes into dissolved organic pool
FracExuC  = 0.7 # Fraction of extra fixed carbon to be exuded

k_sink = 0.05/86400 # Sink rates for particulate carbon and nitrogen (per second)
kDOC   = 1/40/86400 # remineralization rate for DOC, turn over time: 40 days (per second)
kDON   = 1/16/86400 # remineralization rate for DON, turn over time: 16 days (per second)
kPOC   = 1/30/86400 # remineralization rate for POC, turn over time: a month (per second)
kPON   = 1/30/86400 # remineralization rate for PON, turn over time: a month (per second)
####################################### end parameters ##################################
