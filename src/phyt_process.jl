"""
    divide(phyt)
An adult cell divides evenly into two daughter cells
Two daughter cells will be in the same place of the adult cell
"""
function divide(phyt)
    its = size(phyt,1)
    phytos = zeros(Real,its*2)
    # NOT all C and N can turn into new cells
    for i in 1:2
        phytos[1+(i-1)*its]  = phyt[1]          # x
        phytos[2+(i-1)*its]  = phyt[2]          # y
        phytos[3+(i-1)*its]  = phyt[3]          # z
        phytos[4+(i-1)*its]  = phyt[4] .* 0.45  # size
        phytos[5+(i-1)*its]  = phyt[5] .* 0.45  # Bm
        phytos[6+(i-1)*its]  = phyt[6] .* 0.5   # Cq
        phytos[7+(i-1)*its]  = phyt[7] .* 0.5   # Nq
        phytos[8+(i-1)*its]  = phyt[8] .* 0.5   # Pq
        phytos[9+(i-1)*its]  = phyt[9] .* 0.5   # chl
        phytos[10+(i-1)*its] = phyt[10]         # species
        phytos[11+(i-1)*its] = phyt[11].+ 1.0   # generation
        phytos[12+(i-1)*its] = 1.0              # age
        phytos[13+(i-1)*its] = 0.0              # %CfromDOC
        phytos[14+(i-1)*its] = phyt[4] .* 0.45  # age
    end
    return phytos
end

"""
    calc_PAR(surfPAR,grid,Chl,katten_c, katten_w)
Compute PAR for each grid according to depth and Chla concentration
"""
function calc_PAR(surfPAR, grid, Chl, katten_c, katten_w)
    PAR = zeros(grid.Nx, grid.Ny, grid.Nz)
    for k in grid.Nz:-1:1
        atten = zeros(grid.Nx, grid.Ny)
        for i in 1:grid.Nx
            for j in 1:grid.Ny
                atten[i,j] = (Chl[i,j,k] * katten_c + katten_w) * grid.dzF[k]
                PAR[i,j,k] = surfPAR[i,j] * ((1.0 - exp(-atten[i,j])) / atten[i,j])
                surfPAR[i,j] = surfPAR[i,j] * exp(-atten[i,j])
            end
        end
    end
    return PAR
end

"""
    phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
Update the individuals of current time step into next time step
Control flow: Graze ->Grow(photosynthesis, biosynthesis, maintenance) -> Natural Death -> Division
Return a dataframe of next time step individuals, graze number, divide number, death number, and nutrient consumption
"""
function phyt_update(model, ΔT::Int64)
    t = model.t
    clock = t*ΔT%86400÷3600+1 # time of the day, 24-hour
    g = model.grid
    nutrients = model.nutrients
    params = model.params
    phyts_a = copy(model.individuals.phytos)
    idiag = 0

    # load nutrients
    counts = pop_counts(params["P_Nsp"])
    chl_num = count_chl(phyts_a, g)

    # Compute light attenuation, compute from surface
    surfPAR = model.input.PAR[:,:,clock]
    PAR = calc_PAR(surfPAR, g, chl_num, params["katten_c"], params["katten_w"])

    # set up a empty array to record all updated agents
    phyts_b = Real[]
    consume = nutrients_init(g)

    # compute the time index of diagnostics
    diag_t = t*ΔT÷params["diag_freq"]+1

    # iterate phytoplankton agents
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        sp = Int(phyt[10])
        x, y, z = which_grid(phyt, g)
        temp_t = model.input.temp[x,y,z,clock]
        IR_t = PAR[x,y,z]
        NH4 = max(0.0, nutrients.NH4[x, y, z])
        NO3 = max(0.0, nutrients.NO3[x, y, z])
        PO4 = max(0.0, nutrients.PO4[x, y, z])

        # diagnostics
        diag_tmp = zeros(size(params["diag_inds"],1), params["P_Nsp"])
        diag_tmp[1,sp] += 1

        # compute probabilities of grazing
        # Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        if params["Grz_P"] == 0
            P_graz = false
        else
            if (t*ΔT)%3600 ≠ 0 # check hourly
                P_graz = false
            else
                # reg_graz = phyt[4]/params["Grz_P"]
                reg_graz = 1.0/params["Grz_P"]
                P_graz = rand(Bernoulli(reg_graz))
            end
        end

        if P_graz == false #not grazed
            # compute death probability after a certain age, may be abandoned
            reg_age = max(0.0, phyt[12] - params["death_age"])
            shape_factor_death = params["a_death"]*reg_age^params["b_death"]
            P_death = rand(Bernoulli(shape_factor_death/(1+shape_factor_death)))

            if P_death == false # not natural death
                # compute probabilities of division
                P_dvi = false
                if (t*ΔT)%3600 == 0 # check hourly
                    if (params["dvid_type"][sp] == 1) & (phyt[4] ≥ 2.0)
                        reg_size = params["dvid_stp"]*(phyt[4] - params["dvid_size"])
                        reg_divide = 0.2*(tanh(reg_size) + 1)
                        P_dvi = rand(Bernoulli(reg_divide))
                    elseif (params["dvid_type"][sp] == 2) & (phyt[4]-phyt[14]>params["dvid_add"])
                        P_dvi = true
                    end
                end

                if P_dvi == false # not divide
                    # Compute photosynthesis rate
                    α_I = params["α"]*IR_t
                    Tempstd = exp(params["TempAe"]*(1.0/(temp_t+273.15)-1.0/params["Tempref"]))
                    photoTempFunc = params["TempCoeff"]*max(1.0e-10,Tempstd)
                    PCmax_sp = params["PCmax"][sp]/86400
                    PCm = PCmax_sp*photoTempFunc*phyt[4]^params["PC_b"][sp]
                    PC = PCm*(1-exp(-α_I*phyt[9]/(phyt[5]*PCm)))
                    Eₖ = PCm/(phyt[9]/phyt[5]*params["α"])
                    tmp = α_I/params["α"]
                    if (tmp > Eₖ) & (params["inhibcoef"][sp] > 0.0)
                        PC = PC*Eₖ/tmp*params["inhibcoef"][sp]
                    end
                    PS = PC*phyt[5] # unit: mmol C/second/individual
                    PP = PS*ΔT # unit: mmol C/time step/individual
                    # diagnostics
                    diag_tmp[2,sp] += PP

                    # Compute cell-based N uptake rate according Droop limitation
                    Qn = (phyt[7]+phyt[5]*params["R_NC"])/(phyt[5]+phyt[6])
                    #In-Cell N uptake limitation
                    regQn = max(0.0,min(1.0,(params["Nqmax"]-Qn)/(params["Nqmax"]-params["Nqmin"])))
                    VNH4max_sp = params["VNH4max"][sp]/86400
                    VNO3max_sp = params["VNO3max"][sp]/86400
                    VNH4m = VNH4max_sp*phyt[4]^params["VN_b"][sp]
                    VNO3m = VNO3max_sp*phyt[4]^params["VN_b"][sp]
                    NH4uptake = VNH4m*NH4/(NH4+params["KsatNH4"][sp])*regQn
                    NO3uptake = VNO3m*NO3/(NO3+params["KsatNO3"][sp])*regQn
                    VNH4cell = NH4uptake*phyt[5] # unit: mmol N/second/individual
                    VNO3cell = NO3uptake*phyt[5] # unit: mmol N/second/individual
                    VNH4 = min(NH4*g.V[x,y,z]/10.0, VNH4cell*ΔT) # unit: mmol N/time step/individual
                    VNO3 = min(NO3*g.V[x,y,z]/10.0, VNO3cell*ΔT) # unit: mmol N/time step/individual
                    # diagnostics
                    diag_tmp[3,sp] += VNO3
                    diag_tmp[4,sp] += VNH4

                    # Compute cell-based P uptake rate according Droop limitation
                    Qp = (phyt[8]+phyt[5]*params["R_PC"])/(phyt[5]+phyt[6])
                    #In-Cell P uptake limitation
                    regQp = max(0.0,min(1.0,(params["Pqmax"]-Qp)/(params["Pqmax"]-params["Pqmin"])))
                    VPmax_sp = params["VPmax"][sp]/86400
                    VPm = VPmax_sp*phyt[4]^params["VP_b"][sp]
                    Puptake = VPm*PO4/(PO4+params["KsatP"][sp])*regQp
                    VPcell = Puptake*phyt[5] # unit: mmol P/second/individual
                    VPO4 = min(PO4*g.V[x,y,z]/10.0, VPcell*ΔT) # unit: mmol P/time step/individual
                    # diagnostics
                    diag_tmp[5,sp] += VPO4

                    # Compute the ratio of chl synthesis and N uptake
                    # ρ equals to ratio of the realised quantum efficiency for photosynthesis divided by the maximum efficiency
                    if IR_t > 0
                        ρ_chl = PC*params["Chl2N"]/(α_I*phyt[9]/(phyt[5]+phyt[6]))
                    else
                        ρ_chl = 0.0
                    end

                    # Metabolic partitioning for biosynthesis, decrease with size
                    shape_factor_β = params["a_β"]*phyt[4]^params["b_β"]
                    β = shape_factor_β/(1+shape_factor_β)

                    # Compute extra cost for biosynthesis, return a rate (per hour)
                    respir_extra = params["respir_ex"]*phyt[4]^params["respir_b"]

                    # C, N, P storages update
                    phyt[6] = phyt[6] + PP
                    phyt[7] = phyt[7]+ VNH4 + VNO3
                    phyt[8] = phyt[8]+ VPO4

                    # DOC uptake if allowed
                    if params["useDOC"][sp] == 1
                        # read in DOC value at the grid of the individual
                        DOC = max(0.0, nutrients.DOC[x, y, z])
                        # compute the ratio of C reserve to total cellular C
                        Qc = phyt[6] /(phyt[5] + phyt[6])
                        # compute uptake rate of DOC
                        regQc = max(0.0,min(1.0,(params["Cqmax"]-Qc)/(params["Cqmax"]-params["Nqmin"])))
                        VDOCmax_sp = params["VDOCmax"][sp]/86400
                        VDOCm = VDOCmax_sp*phyt[4]^params["VDOC_b"][sp]
                        DOCuptake = VDOCm*DOC/(DOC+params["KsatDOC"][sp])*regQc
                        VDOCcell = DOCuptake*phyt[5] # unit: mmol C/second/individual
                        VDOC = min(DOC*g.V[x,y,z]/10.0, VDOCcell*ΔT) # unit: mmol C/time step/individual
                        # diagnostics
                        diag_tmp[6,sp] += VDOC
                        # update C reserve of the individual
                        phyt[6] = phyt[6] + VDOC
                        # add up consume of DOC by DOC uptake
                        consume.DOC[x, y, z] = consume.DOC[x, y, z] - VDOC
                        # compute percentage of C that is from DOC uptake
                        phyt[13] = VDOC/(VDOC+PP)
                    end

                    # maximum biosynthesis rate based on carbon availability
                    k_mtb = params["k_mtb"]*phyt[4]^params["b_k_mtb"]
                    BS_Cmax = β*k_mtb*ΔT*phyt[6]/(1+respir_extra)
                    MaintenC = (1-β)*BS_Cmax/β

                    # maximum allowed biosynthesis rate by Nq and Pq
                    BS_Nmax = k_mtb*ΔT*phyt[7]/params["R_NC"]
                    BS_Pmax = k_mtb*ΔT*phyt[8]/params["R_PC"]

                    # acutall biosynthesis rate & excretion
                    BS_C = min(BS_Cmax, BS_Nmax, BS_Pmax)
                    excretC = max(0.0, BS_Cmax-BS_C)

                    # diagnostics
                    diag_tmp[7,sp] += BS_C
                    diag_tmp[8,sp] += MaintenC
                    diag_tmp[9,sp] += excretC

                    # update quotas, biomass, Chla and cell size etc.
                    CostC = BS_C*(1+respir_extra)
                    phyt[5] = phyt[5] + BS_C
                    phyt[6] = phyt[6] - CostC -MaintenC - excretC
                    phyt[7] = phyt[7]- BS_C*params["R_NC"]
                    phyt[8] = phyt[8]- BS_C*params["R_PC"]
                    dsize= BS_C/(params["P_Cquota"][sp]*params["P_Nsuper"]) # normalized by standard C quota
                    phyt[4] = max(0.0,phyt[4]+dsize)
                    phyt[9] = phyt[9] + ρ_chl*BS_C*params["R_NC"]
                    phyt[12]= phyt[12] + 1.0*(ΔT/3600)

                    # diagnostics
                    diag_tmp[13,sp] += phyt[5]
                    diag_tmp[14,sp] += phyt[6]
                    diag_tmp[15,sp] += phyt[7]
                    diag_tmp[16,sp] += phyt[8]
                    diag_tmp[17,sp] += phyt[9]

                    consume.DIC[x, y, z] = consume.DIC[x, y, z] + MaintenC + CostC - BS_C
                    consume.DOC[x, y, z] = consume.DOC[x, y, z] + excretC
                    consume.NH4[x, y, z] = consume.NH4[x, y, z] - VNH4
                    consume.NO3[x, y, z] = consume.NO3[x, y, z] - VNO3
                    consume.PO4[x, y, z] = consume.PO4[x, y, z] - VPO4
                    append!(phyts_b,phyt)
                else # divide
                    counts.divid[sp] += 1
                    diag_tmp[10,sp] += 1
                    phyts = divide(phyt)
                    append!(phyts_b,phyts)
                    consume.DIC[x, y, z] = consume.DIC[x, y, z] + phyt[5]*0.1 # consume C when cell is divided
                end # divide
            else # natural death
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt[5]+phyt[6])*params["mortFracC"]
                consume.DON[x, y, z] = consume.DON[x, y, z] + (phyt[7]+phyt[5]*params["R_NC"])*params["mortFracN"]
                consume.DOP[x, y, z] = consume.DOP[x, y, z] + (phyt[8]+phyt[5]*params["R_PC"])*params["mortFracP"]
                consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt[5]+phyt[6])*(1.0 - params["mortFracC"])
                consume.PON[x, y, z] = consume.PON[x, y, z] + (phyt[7]+phyt[5]*params["R_NC"])*(1.0 - params["mortFracN"])
                consume.POP[x, y, z] = consume.POP[x, y, z] + (phyt[8]+phyt[5]*params["R_PC"])*(1.0 - params["mortFracP"])
                counts.death[sp] += 1
                diag_tmp[12,sp] += 1
            end # naturan death
        else #grazed, no sloppy feeding here, all nutrients go back to organic pools
            counts.graze[sp] += 1
            diag_tmp[11,sp] += 1
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt[5]+phyt[6])*params["grazFracC"]
            consume.DON[x, y, z] = consume.DON[x, y, z] + (phyt[7]+phyt[5]*params["R_NC"])*params["grazFracN"]
            consume.DOP[x, y, z] = consume.DOP[x, y, z] + (phyt[8]+phyt[5]*params["R_PC"])*params["grazFracP"]
            consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt[6]+phyt[6])*(1.0 - params["grazFracC"])
            consume.PON[x, y, z] = consume.PON[x, y, z] + (phyt[7]+phyt[5]*params["R_NC"])*(1.0 - params["grazFracN"])
            consume.POP[x, y, z] = consume.POP[x, y, z] + (phyt[8]+phyt[5]*params["R_PC"])*(1.0 - params["grazFracP"])
        end # graze
        # diagnostics
        for it in size(params["diag_inds"],1)
            if params["diag_inds"][it] == 1
                idiag += 1
                model.diags.spcs[x,y,z,diat_t,sp,idiag] += diag_tmp[it,sp]
            end
        end
    end # for loop to traverse the array of agents
    # diagnostics
    model.diags.tr[:,:,:,diag_t,1] = PAR
    phyts_b = reshape(phyts_b,size(phyts_a,1),Int(length(phyts_b)/size(phyts_a,1)))
    return phyts_b,counts,consume
end # for loop of time
