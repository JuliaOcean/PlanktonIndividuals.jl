var documenterSearchIndex = {"docs":
[{"location":"equations/#parameterization-of-phytoplankton-physiology-1","page":"Equations","title":"parameterization of phytoplankton physiology","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"(Image: skematic)","category":"page"},{"location":"equations/#Photosynthesis-1","page":"Equations","title":"Photosynthesis","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"Basically, we follow Geider et al. 1998 for the parameterization of photosynthesis, but without nutrient limitation.","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"PP=PP_maxcdot (1-e^frac-alpha cdot Icdot ChlPP_maxcdot C)","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"Unit: mmolC/cell/s","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"PP_max is scaled by a power-law relationship of cell size","category":"page"},{"location":"equations/#Nutrient-Uptake-1","page":"Equations","title":"Nutrient Uptake","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"Include intracellular nutrient limitation:","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"RegQ_i=biggfracR_iC^max-Q_iR_iC^max-R_iC^minbigg_0^1\nV_i=V_i^maxcdot regQ_icdotfracii+K_i^sat","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"i denotes NH_4, NO_3, PO_4.","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"Unit: mmolN/cell/s","category":"page"},{"location":"equations/#Biosynthesis-and-Excretion-1","page":"Equations","title":"Biosynthesis & Excretion","text":"","category":"section"},{"location":"equations/#Update-reserves:-1","page":"Equations","title":"Update reserves:","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"Q_C^R = Q_C^R+PP\nQ_N^R=Q_N^R+V_NO_3+V_NH_4\nQ_P^R=Q_P^R+V_PO_4","category":"page"},{"location":"equations/#Metabolic-Partitioning-1","page":"Equations","title":"Metabolic Partitioning","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"beta=fracacdot Vol_cell^b1+acdot Vol_cell^b\nBioSynC = betacdot k_mtbcdot Q_C^R\nMaintenC=(1-beta)cdot k_mtbcdot Q_C^R","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"BioSynN = k_mtbcdot Q_N^RR_NC\nBioSynP = k_mtbcdot Q_P^RR_PC","category":"page"},{"location":"equations/#Compute-biosynthesis-rate-and-excretion-rate-1","page":"Equations","title":"Compute biosynthesis rate and excretion rate","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"BioSyn=min(BioSynCBioSynNBioSynP)\nExcretC=BioSynC - BioSyn","category":"page"},{"location":"equations/#Update-reserves-and-biomass-1","page":"Equations","title":"Update reserves and biomass","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"Q_C^B = Q_C^B + BioSyn\nQ_C^R = Q_C^R - BioSyn - MaintenC\nQ_N^R = Q_N^R - BioSyn*R_NC\nQ_P^R = Q_N^R - BioSyn*R_PC","category":"page"},{"location":"equations/#Cell-division-1","page":"Equations","title":"Cell division","text":"","category":"section"},{"location":"equations/#","page":"Equations","title":"Equations","text":"We use relative cell size RCS to indicate cell division. RCS of the smallest cell is 1.0. Q_C^B of the smallest cell is Cquota. Cells start to divide at RCS=20. The probability of individual cell division is a sigmoidal function of RCS.","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"RCS = Q_C^B  Cquota\nP_divide = rand(Bernoulli(02*(tanh(a*(RCS-b))+1)))","category":"page"},{"location":"equations/#","page":"Equations","title":"Equations","text":"P_divide is computed every hour no matter what time step the model is.","category":"page"},{"location":"#PlanktonIndividuals.jl-1","page":"Home","title":"PlanktonIndividuals.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Documentation for PlanktonIndividuals.jl which simulates the behavior of an ensemble of phytoplankton individuals.","category":"page"},{"location":"#Use-Example-1","page":"Home","title":"Use Example","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Here we use Oceananigans.jl to generate velocity fields and then use those to drive the individual-based model.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pkg.develop(PackageSpec(path=\"PlanktonIndividuals.jl\"))\nusing PlanktonIndividuals\np = dirname(pathof(PlanktonIndividuals))\ninclude(joinpath(p,\"Oceananigans_PlanktonIndividuals.jl\"))","category":"page"},{"location":"#Unit-Testing-1","page":"Home","title":"Unit Testing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"The tests use input files from samples/. The test suite includes zero-, one-, two-, and three-dimensional simulations.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pkg.develop(PackageSpec(path=\"PlanktonIndividuals.jl\"))\nPkg.test(\"PlanktonIndividuals\")","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Contents:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"index.md\",\"various.md\"]\nDepth = 3","category":"page"},{"location":"#API-Guide-1","page":"Home","title":"API Guide","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [PlanktonIndividuals]\nOrder   = [:type,:function]","category":"page"},{"location":"#PlanktonIndividuals.PI_Model-Tuple{Architecture,Any,Any}","page":"Home","title":"PlanktonIndividuals.PI_Model","text":"PI_model(arch, grid, RunParam; t, nutrients, individuals, PARF, temp, params, diag_ntrs, diag_nprocs)\n\nGenerate the PI_Model struct on grid. \n\nKeyword Arguments\n\narch (required): CPUs() or GPUs(). The computer architecture used to time-step model.\ngrid (required): The resolution and discrete geometry on which model is solved.\nRunParam (required): Run time parameters for model including number of time steps, ΔT, etc.\nnutrients (required): Nutrient tracers initialized with grid and initial conditions.\nindividuals (optional):  Individuals generated according to RunParam.\nPARF and temp (optional): External forcings of PAR and temperature.\nparams (optional): Parameter set for biogeochemical processes modeled in the model.\ndiag_ntrs (optional): NamedTuple containing the names of nutrient fields to be diagnosed\ndiag_nprocs (optional): NamedTuple containing the names of physiological processes to be diagnosed\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.PI_TimeStep!-Tuple{PI_Model,Any,String}","page":"Home","title":"PlanktonIndividuals.PI_TimeStep!","text":"PI_TimeStep!(model, ΔT, resultpath)\n\nUpdate physiology processes and nutrient field of PI_Model one time step forward.\n\nKeyword Arguments\n\nmodel: PI_Model to be updated one time step forward\nΔT: The length of a time step\nresultpath (optional): The file path to store model output. \n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.PrepRunDir","page":"Home","title":"PlanktonIndividuals.PrepRunDir","text":"PrepRunDir(res::String=\"results/\")\n\nCreate res/ folder if needed. Remove old files from it if needed.\n\n\n\n\n\n","category":"function"},{"location":"#PlanktonIndividuals.gen_Grid-Tuple{}","page":"Home","title":"PlanktonIndividuals.gen_Grid","text":"gen_Grid(size = (Nx, Ny, Nz), spacing = (Δx, Δy, Δz), halo = (2, 2, 2))\n\nCreats a Grids struct with size = (Nx, Ny, Nz) grid points.\n\nKeyword Arguments\n\nsize (required): A tuple prescribing the number of grid points.                        size is a 3-tuple no matter for 3D, 2D, or 1D model.\nspacing (required): A tuple prescribing the length of each grid point in x, y, and z directions.                       spacing is a 3-tuple no matter for 3D, 2D, or 1D model.\nhalo (optional): A tuple of integers that specifies the size of the halo region of cells                       surrounding the physical interior for each direction.                       halo is a 3-tuple no matter for 3D, 2D, or 1D model.\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.gen_nutrients-Tuple{Any,Any,Any}","page":"Home","title":"PlanktonIndividuals.gen_nutrients","text":"gen_nutrients(arch, grid, nut)\n\nSet up initial nutrient fields according to grid.\n\nKeyword Arguments\n\narch: CPUs() or GPUs(). The computer architecture used to time-step model.\ngrid: The resolution and discrete geometry on which nutrient fields are solved.\nnut: An 10-element array with each element representing the initial condition of a kind of nutrient.\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.load_nut_initials-Tuple{Any,Any,Any}","page":"Home","title":"PlanktonIndividuals.load_nut_initials","text":"load_nut_initials(arch, paths, grid)\n\nLoad nutrient initial conditions from files\n\nKeyword Arguments\n\narch: CPUs() or GPUs(). The computer architecture used to time-step model.\npaths: NamedTuple containing the file paths pointing to the files of nutrient initial conditions.\ngrid: The resolution and discrete geometry on which nutrient fields are solved.\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.update_params!-Tuple{Dict,Dict}","page":"Home","title":"PlanktonIndividuals.update_params!","text":"update_params!(parameters, tmp)\n\nUpdate parameter values based on a .yaml file provided by user\n\nKeyword Arguments\n\nparameters is default parameter set\ntmp is the parameters read from .yaml file and needs to update\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.write_diags_to_jld2-NTuple{4,Any}","page":"Home","title":"PlanktonIndividuals.write_diags_to_jld2","text":"write_diags_to_jld2(diags, filepath, t, ncounts)\n\nwrite model output of individuals at each time step to a binary file\n\nKeyword Arguments\n\ndiags: NamedTuple of a list of diagnostics at current time step.\nfilepath: The file path to store JLD2 files.\nt: Current time of model in second, usually starting from 0.\nncounts: the number of time steps included in each diagnostic\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.write_individuals_to_jld2-Tuple{NamedTuple,Any,Any}","page":"Home","title":"PlanktonIndividuals.write_individuals_to_jld2","text":"write_individuals_to_bin(phytos, filepath, t)\n\nwrite model output of individuals at each time step to a binary file\n\nKeyword Arguments\n\nphytos: NamedTuple of a list of individual species.\nfilepath: The file path to store JLD2 files.\nt: Current time of model in second, usually starting from 0.\natts (optional): attributes of individuals to save, default (:x, :y, :z)\n\n\n\n\n\n","category":"method"},{"location":"#PlanktonIndividuals.write_nut_nc_each_step-Tuple{NamedTuple,Int64,String}","page":"Home","title":"PlanktonIndividuals.write_nut_nc_each_step","text":"write_nut_nc_each_step(nut, t, filepath)\n\nWrite a NetCDF file of nutrient fields at each time step\n\nKeyword Arguments\n\nnut: NamedTuple of nutrient tracers.\nt: Current time of model in second, usually starting from 0.\nfilepath: The file path to store NetCDF files.\n\n\n\n\n\n","category":"method"},{"location":"various/#Function-Inventory-1","page":"Various","title":"Function Inventory","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"Use cases:","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"Oceananigans_PlanktonIndividuals.jl","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"Infrastructure functions:","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"model_struct.jl defines velocity, grids, nutrient_fields, and rem structs.\nmodel.jl deals with the model basic functionalities.\nmodel_setup.jl is setup_agents and setup_nutrients\noption_params.jl and param_defaults.jl deal with model parameter values\nutils.jl utility functions for IO and operations.","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"process functions:","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"phyt_process.jl = daynight, PAR_cal, PC, Nuptake, chl_sync, divide, phyt_update\nnutrient_processes.jl = compute_nut_biochem, compute_source_term, nut_update\ndst3fl.jl 3rd order DST Scheme with flux limiting\n2nd_adv_diffu.jl right hand side term functions (?)\nagent_div.jl = is mostly agent_move, agent_move_1D + double_grid, trilinear_itpl, simple_itpl","category":"page"},{"location":"various/#Model-Variables-1","page":"Various","title":"Model Variables","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"The lists and text below likely is out of date","category":"page"},{"location":"various/#)-state-variables-1","page":"Various","title":"1) state variables","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"phyts_a\nnutrients","category":"page"},{"location":"various/#)-input-variables-1","page":"Various","title":"2) input variables","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"g from grid_offline() (grid variables)\ntemp, IR from read_input() (before loop)\nremin from rem() (before loop)\nvelᵇ from read_offline_vels() (in loop)","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"Notes:","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"IR and temp get cycle through / interpolated inside phyt_update (via IR[trunc(Int,t*ΔT/3600)]). Not sure about daynight(t,IR) in phyt_update. \nread_offline_vels receives trunc(Int,t*ΔT/3600) as argument.\nremin is time-invariant; like g.","category":"page"},{"location":"various/#)-intermediate-variables-1","page":"Various","title":"3) intermediate variables","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"F = compute_nut_biochem(nutrients, remin)\ngtr = compute_source_term(nutrients, velᵇ, g, F)\nnutₜ = nut_update(nutrients, consume, g, gtr, ΔT)","category":"page"},{"location":"various/#)-diagnostic-variables-1","page":"Various","title":"4) diagnostic variables","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"dvid_ct, graz_ct, death_ct from phyt_update\ngtr, nutₜ, velᵇ, agent_num from ","category":"page"},{"location":"various/#Output-Files-1","page":"Various","title":"Output Files","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"The lists and text below likely is out of date","category":"page"},{"location":"various/#)-listing-1","page":"Various","title":"1) listing","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"B1.bin\t\tall agents at all time steps for species 1\nB2.bin\t\t... species 2\ncons_C.txt\ttotal carbon for each time step\ncons_DIN.txt\t... DIN ...\ncons_N.txt\t... N ...\ngrid.bin\tgrid (input)\nIR.bin\t\tirradiance (input)\nnutrients/nut.0001.nc (all nitrients for 1 time step)\t\nnutrients/nut.0002.nc (same...)\noutput1.bin\taverage for all agents for each time step for species 1\noutput2.bin\t... 2\noutput.bin\t... for the two species\nVD1.bin\t\tvertical profile of the agents opouplation for species 1\nVD2.bin\t\t... species 2","category":"page"},{"location":"various/#)-netcdf-output-1","page":"Various","title":"2) netcdf output","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"netcdf nut.0001 {\ndimensions:\n        yC = 1 ;\n        xC = 1 ;\n        zC = 40 ;\nvariables:\n        double DIC(zC, yC, xC) ;\n                DIC:units = \"mmolC/m^3\" ;\n        float yC(yC) ;\n                yC:units = \"m\" ;\n                yC:longname = \"Locations of the cell centers in the y-direction.\" ;\n        float xC(xC) ;\n                xC:units = \"m\" ;\n                xC:longname = \"Locations of the cell centers in the x-direction.\" ;\n        float zC(zC) ;\n                zC:units = \"m\" ;\n                zC:longname = \"Locations of the cell centers in the z-direction.\" ;\n        double DIN(zC, yC, xC) ;\n                DIN:units = \"mmolN/m^3\" ;\n        double DOC(zC, yC, xC) ;\n                DOC:units = \"mmolC/m^3\" ;\n        double DON(zC, yC, xC) ;\n                DON:units = \"mmolN/m^3\" ;\n        double POC(zC, yC, xC) ;\n                POC:units = \"mmolC/m^3\" ;\n        double PON(zC, yC, xC) ;\n                PON:units = \"mmolN/m^3\" ;","category":"page"},{"location":"various/#Data-Structures-1","page":"Various","title":"Data Structures","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"The lists and text below likely is out of date","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":" output = DataFrame(time=0, \n gen_ave=mean(B[1].gen), \n spec_ave = mean(B[1].sp),\n Cq1_ave=mean(B[1].Cq1), \n Cq2_ave=mean(B[1].Cq2), \n Nq_ave=mean(B[1].Nq),\n size_ave=mean(B[1].size),\n chl_ave=mean(B[1].chl),\n Population=size(B[1],1),\n dvid=0,\n graz=0,\n death = 0)","category":"page"},{"location":"various/#Plotting-Functions-1","page":"Various","title":"Plotting Functions","text":"","category":"section"},{"location":"various/#","page":"Various","title":"Various","text":"The lists and text below likely is out of date","category":"page"},{"location":"various/#","page":"Various","title":"Various","text":"julia> using DataFrames, Serialization\njulia> output=open(deserialize,\"results/output.bin\");\njulia> plot(output.Population)","category":"page"}]
}
