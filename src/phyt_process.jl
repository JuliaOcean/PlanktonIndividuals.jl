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
        phytos[13+(i-1)*its] = phyt[4] .* 0.45  # init size
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
    clock = t%86400÷ΔT+1 # time of the day, in number of time steps
    g = model.grid
    nutrients = model.nutrients
    params = model.params
    phyts_a = copy(model.individuals.phytos)

    # load nutrients
    chl_num = count_chl(phyts_a, g)

    # Compute light attenuation, compute from surface
    surfPAR = model.input.PAR[:,:,clock]
    PAR = calc_PAR(surfPAR, g, chl_num, params["katten_c"], params["katten_w"])

    # set up a empty array to record all updated agents
    phyts_b = Real[]
    consume = nutrients_init(g)

    # compute the time index of diagnostics
    diag_t = t÷params["diag_freq"]+1
    npop = zeros(Int, g.Nx, g.Ny, g.Nz, params["P_Nsp"])

    num_phyt = size(phyts_a,2)

    # iterate phytoplankton agents
    for i in 1:num_phyt
        phyt = phyts_a[:,i]
        sp = Int(phyt[10])
        x, y, z = which_grid(phyt, g)
        temp_t = model.input.temp[x,y,z,clock]
        IR_t = PAR[x,y,z]
        NH4 = max(0.0, nutrients.NH4[x, y, z])
        NO3 = max(0.0, nutrients.NO3[x, y, z])
        PO4 = max(0.0, nutrients.PO4[x, y, z])

        #diagnostics
        idiag = 0
        npop[x, y, z, sp] += 1

        # compute probabilities of grazing
        # Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        if params["Grz_P"] == 0
            P_graz = false
        else
            if t%600 ≠ 1 # check every 10 mins
                P_graz = false
            else
                reg_graz = 1.0/params["Grz_P"]
                # reg_graz = exp(num_phyt/(params["P_Nind"]*params["P_Nsp"]))/params["Grz_P"]
                P_graz = rand(Bernoulli(reg_graz))
            end
        end

        if P_graz == false #not grazed
            # compute death probability according to cell size
            # minimal vital cell size is 1.0
            P_death = false
            if t%600 == 1 # check every 10 mins
                reg_de = 6.0*(params["death_reg"][sp] - phyt[4])
                reg_death = params["P_death"][sp]*(tanh(reg_de) + 1)
                P_death = rand(Bernoulli(reg_death))
            end

            if P_death == false # not natural death
                # compute probabilities of division
                P_dvi = false
                if t%600 == 1 # check every 10 mins
                    if params["dvid_type"][sp] == 1 # sizer-like cell division
                        if phyt[5] ≥ 2*params["P_Cquota"][sp]*params["P_Nsuper"]
                            reg_size   = params["dvid_stp"][sp]*(phyt[4] - params["dvid_size"][sp])
                            reg_divide = params["P_dvid"][sp]*(tanh(reg_size) + 1)
                            P_dvi      = rand(Bernoulli(reg_divide))
                        end
                    elseif params["dvid_type"][sp] == 2 # adder-like cell division
                        add_size = phyt[4] - phyt[13]
                        if phyt[5] ≥ 2*params["P_Cquota"][sp]*params["P_Nsuper"]
                            reg_add    = params["dvid_stp"][sp]*(add_size - params["dvid_add"][sp])
                            reg_divide = params["P_dvid"][sp]*(tanh(reg_add) + 1)
                            P_dvi      = rand(Bernoulli(reg_divide))
                        end
                    elseif params["dvid_type"][sp] == 3 # timer-like (age) cell division
                        if phyt[5] ≥ 2*params["P_Cquota"][sp]*params["P_Nsuper"]
                            reg_age    = params["dvid_stp"][sp]*(phyt[12] - params["dvid_age"][sp])
                            reg_divide = params["P_dvid"][sp]*(tanh(reg_age) + 1)
                            P_dvi      = rand(Bernoulli(reg_divide))
                        end
                    elseif params["dvid_type"][sp] == 4 # timer-like (circadian clock) cell division
                        if phyt[5] ≥ 2*params["P_Cquota"][sp]*params["P_Nsuper"]
                            cirT       = t % 86400 ÷ 3600
                            reg_cirT   = params["dvid_stp"][sp]*(cirT - params["dvid_cirT"][sp])
                            reg_divide = params["P_dvid"][sp]*(tanh(reg_cirT) + 1)
                            P_dvi      = rand(Bernoulli(reg_divide))
                        end
                    elseif params["dvid_type"][sp] == 5 # timer-like (circadian clock) cell division
                        if phyt[5] ≥ 2*params["P_Cquota"][sp]*params["P_Nsuper"]
                            # use light intensity to indicate circadian clock in the cell
                            reg_par    = params["dvid_stp"][sp]*(params["dvid_par"][sp] - IR_t)
                            reg_divide = params["P_dvid"][sp]*(tanh(reg_par) + 1)
                            P_dvi      = rand(Bernoulli(reg_divide))
                        end
                    elseif params["dvid_type"][sp] == 6 # timer & sizer-like cell division
                        cirT       = t % 86400 ÷ 3600
                        reg_cirT   = cirT - params["dvid_cirT"][sp]
                        reg_size   = params["dvid_stp"][sp]*(phyt[4] - params["dvid_size"][sp])
                        reg_divide = params["P_dvid"][sp]*(tanh(reg_size) + 1)*(tanh(reg_cirT) + 1)/2.0
                        P_dvi      = rand(Bernoulli(reg_divide))
                    elseif params["dvid_type"][sp] == 7 # timer-like (circadian clock) cell division
                        cirT       = t % 86400 ÷ 3600
                        reg_cirT   = params["dvid_stp"][sp]*(cirT - params["dvid_cirT"][sp])
                        reg_divide = params["P_dvid"][sp]*(tanh(reg_cirT) + 1)
                        P_dvi      = rand(Bernoulli(reg_divide))
                    elseif params["dvid_type"][sp] == 8 # sizer-like cell division
                        reg_size   = params["dvid_stp"][sp]*(phyt[4] - params["dvid_size"][sp])
                        reg_divide = params["P_dvid"][sp]*(tanh(reg_size) + 1)
                        P_dvi      = rand(Bernoulli(reg_divide))
                    else
                        print("wrong division type! \n")
                    end
                end

                if P_dvi == false # not divide
                    # diagnostics
                    if params["diag_inds"][1] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += 1
                    end
                    # Compute photosynthesis rate
                    α_I = params["α"][sp]*IR_t*params["Φ"][sp]
                    Tempstd = exp(params["TempAe"]*(1.0/(temp_t+273.15)-1.0/params["Tempref"]))
                    photoTempFunc = params["TempCoeff"]*max(1.0e-10,Tempstd)
                    PCmax_sp = params["PCmax"][sp]/86400
                    PCm = PCmax_sp*photoTempFunc*phyt[4]^params["PC_b"][sp]
                    PC = PCm*(1-exp(-α_I*phyt[9]/(phyt[5]*PCm)))
                    # photo-inhibition
                    if params["inhibcoef"][sp] > 0.0
                        Eₖ = PCm/(phyt[9]/phyt[5]*params["α"][sp]*params["Φ"][sp])
                        tmp = α_I/(params["α"][sp]*params["Φ"][sp])
                        if tmp > Eₖ
                            PC = PC*Eₖ/tmp*params["inhibcoef"][sp]
                        end
                    end
                    PS = PC*phyt[5] # unit: mmol C/second/individual
                    PP = PS*ΔT # unit: mmol C/time step/individual
                    #diagnostics
                    if params["diag_inds"][2] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += PP
                    end

                    if params["isDiaz"] ≠ 1
                        # Compute cell-based N uptake rate according Droop limitation
                        Qn = (phyt[7]+phyt[5]*params["R_NC"])/(phyt[5]+phyt[6])
                        #In-Cell N uptake limitation
                        regQn = max(0.0,min(1.0,(params["Nqmax"][sp]-Qn)/(params["Nqmax"][sp]-params["Nqmin"][sp])))
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
                    else
                        VNO3 = 0.0
                        VNH4max_sp = params["VNH4max"][sp]/86400
                        VNH4 = VNH4max_sp*phyt[4]^params["VN_b"][sp]*phyt[5]*ΔT # unit: mmol N/time step/individual
                        # cost energy
                        phyt[6] = phyt[6] - VNH4*params["C2Nfix"][sp]
                        if phyt[6] ≤ 0.0
                            exceed = 0.0 - phyt[6]
                            VNH4 = VNH4 - exceed/params["C2Nfix"][sp]
                            phyt[6] = 0.0
                        end
                    end
                    #diagnostics
                    if params["diag_inds"][3] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += VNO3
                    end
                    if params["diag_inds"][4] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += VNH4
                    end

                    # Compute cell-based P uptake rate according Droop limitation
                    Qp = (phyt[8]+phyt[5]*params["R_PC"])/(phyt[5]+phyt[6])
                    #In-Cell P uptake limitation
                    regQp = max(0.0,min(1.0,(params["Pqmax"][sp]-Qp)/(params["Pqmax"][sp]-params["Pqmin"][sp])))
                    VPmax_sp = params["VPmax"][sp]/86400
                    VPm = VPmax_sp*phyt[4]^params["VP_b"][sp]
                    Puptake = VPm*PO4/(PO4+params["KsatP"][sp])*regQp
                    VPcell = Puptake*phyt[5] # unit: mmol P/second/individual
                    VPO4 = min(PO4*g.V[x,y,z]/10.0, VPcell*ΔT) # unit: mmol P/time step/individual
                    #diagnostics
                    if params["diag_inds"][5] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += VPO4
                    end

                    # Compute the ratio of chl synthesis and N uptake
                    # ρ equals to ratio of the realised quantum efficiency for photosynthesis divided by the maximum efficiency
                    if IR_t > 0
                        ρ_chl = PC*params["Chl2N"]/(α_I*phyt[9]/(phyt[5]))
                    else
                        ρ_chl = 0.0
                    end

                    # # Metabolic partitioning for biosynthesis, decrease with size
                    # shape_factor_β = params["a_β"][sp]*phyt[4]^params["b_β"][sp]
                    # β = shape_factor_β/(1+shape_factor_β)

                    # Respiration
                    respir_sp = params["respir_a"][sp]/86400
                    respir_m = respir_sp*phyt[4]^params["respir_b"][sp]*photoTempFunc
                    respir_cell = respir_m*phyt[5] # unit: mmolC/second/individual
                    respir = respir_cell*ΔT # unit: mmolC/time step/individual

                    # DOC uptake if allowed
                    if params["useDOC"][sp] == 1
                        # read in DOC value at the grid of the individual
                        DOC = max(0.0, nutrients.DOC[x, y, z])
                        # compute the ratio of C reserve to total cellular C
                        Qc = phyt[6] /(phyt[5] + phyt[6])
                        # compute uptake rate of DOC
                        regQc = max(0.0,min(1.0,(params["Cqmax"][sp]-Qc)/(params["Cqmax"][sp]-params["Cqmin"][sp])))
                        VDOCmax_sp = params["VDOCmax"][sp]/86400
                        VDOCm = VDOCmax_sp*phyt[4]^params["VDOC_b"][sp]
                        DOCuptake = VDOCm*DOC/(DOC+params["KsatDOC"][sp])*regQc
                        VDOCcell = DOCuptake*phyt[5] # unit: mmol C/second/individual
                        VDOC = min(DOC*g.V[x,y,z]/10.0, VDOCcell*ΔT) # unit: mmol C/time step/individual
                    else
                        VDOC = 0.0
                    end
                    #diagnostics
                    if params["diag_inds"][6] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += VDOC
                    end

                    # C, N, P storages update
                    phyt[6] = phyt[6] + PP + VDOC - respir
                    if phyt[6] ≤ 0.0 # if C reserve is not enough for respiration
                        exceed = 0.0 - phyt[6]
                        phyt[5] = phyt[5] - exceed # respire functional C instead
                        phyt[6] = 0.0
                        phyt[7] = phyt[7] + exceed*params["R_NC"] # return N from function pool to N reserve
                        phyt[8] = phyt[8] + exceed*params["R_PC"] # return P from function pool to P reserve
                    end
                    phyt[7] = phyt[7] + VNH4 + VNO3
                    phyt[8] = phyt[8] + VPO4

                    # maximum biosynthesis rate based on carbon availability
                    k_mtb = params["k_mtb"][sp]*phyt[4]^params["b_k_mtb"][sp]/86400 # per second
                    # BS_Cmax = β*k_mtb*ΔT*phyt[6]/(1+respir_extra)
                    # MaintenC = (1-β)*BS_Cmax/β
                    BS_Cmax = k_mtb*phyt[6]*ΔT

                    # maximum allowed biosynthesis rate by Nq and Pq
                    BS_Nmax = k_mtb*ΔT*phyt[7]/params["R_NC"]
                    BS_Pmax = k_mtb*ΔT*phyt[8]/params["R_PC"]

                    # acutall biosynthesis rate & excretion
                    BS_C = min(BS_Cmax, BS_Nmax, BS_Pmax)
                    excretC = max(0.0, BS_Cmax-BS_C)

                    #diagnostics
                    if params["diag_inds"][7] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += BS_C
                    end
                    if params["diag_inds"][8] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += respir
                    end
                    if params["diag_inds"][9] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += excretC
                    end

                    # update quotas, biomass, Chla and cell size etc.
                    phyt[5] = phyt[5] + BS_C
                    phyt[6] = phyt[6] - BS_C - excretC
                    phyt[7] = phyt[7] - BS_C*params["R_NC"]
                    phyt[8] = phyt[8] - BS_C*params["R_PC"]
                    phyt[9] = phyt[9] + ρ_chl*BS_C*params["R_NC"]
                    phyt[12]= phyt[12]+ 1.0*(ΔT/3600)

                    # normalized by standard C quota
                    # dsize= (PP + VDOC - MaintenC - excretC)/(params["P_Cquota"][sp]*params["P_Nsuper"])
                    # phyt[4] = max(0.0,phyt[4]+dsize)
                    phyt[4] = (phyt[5] + phyt[6])/(params["P_Cquota"][sp]*params["P_Nsuper"])

                    #diagnostics
                    if params["diag_inds"][10] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += phyt[5]
                    end
                    if params["diag_inds"][11] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += phyt[6]
                    end
                    if params["diag_inds"][12] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += phyt[7]
                    end
                    if params["diag_inds"][13] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += phyt[8]
                    end
                    if params["diag_inds"][14] == 1
                        idiag += 1
                        model.diags.spcs[x,y,z,diag_t,sp,idiag] += phyt[9]
                    end

                    consume.DIC[x, y, z] = consume.DIC[x, y, z] + respir + BS_C - PP
                    consume.DOC[x, y, z] = consume.DOC[x, y, z] + excretC - VDOC
                    consume.NH4[x, y, z] = consume.NH4[x, y, z] - VNH4
                    consume.NO3[x, y, z] = consume.NO3[x, y, z] - VNO3
                    consume.PO4[x, y, z] = consume.PO4[x, y, z] - VPO4
                    append!(phyts_b,phyt)
                else # divide
                    model.diags.pop[x,y,z,diag_t,sp,1] += 1
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
                model.diags.pop[x,y,z,diag_t,sp,2] += 1
            end # naturan death
        else #grazed, no sloppy feeding here, all nutrients go back to organic pools
            model.diags.pop[x,y,z,diag_t,sp,3] += 1
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt[5]+phyt[6])*params["grazFracC"]
            consume.DON[x, y, z] = consume.DON[x, y, z] + (phyt[5]*params["R_NC"]+phyt[7])*params["grazFracN"]
            consume.DOP[x, y, z] = consume.DOP[x, y, z] + (phyt[5]*params["R_PC"]+phyt[8])*params["grazFracP"]
            consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt[5]+phyt[6])*(1.0 - params["grazFracC"])
            consume.PON[x, y, z] = consume.PON[x, y, z] + (phyt[5]*params["R_NC"]+phyt[7])*(1.0 - params["grazFracN"])
            consume.POP[x, y, z] = consume.POP[x, y, z] + (phyt[5]*params["R_PC"]+phyt[8])*(1.0 - params["grazFracP"])
        end # graze
    end # for loop to traverse the array of agents
    # diagnostics
    model.diags.tr[:,:,:,diag_t,1] += PAR
    for k in 1:params["P_Nsp"]
        model.diags.tr[:,:,:,diag_t,k+1] += npop[:,:,:,k]
    end
    phyts_b = reshape(phyts_b,size(phyts_a,1),Int(length(phyts_b)/size(phyts_a,1)))
    return phyts_b,consume
end # for loop of time
