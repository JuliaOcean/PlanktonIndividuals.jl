###########################################
# define functions of pysiology processes #
###########################################
"""
    divide(phyt)
An adult cell divides evenly into two daughter cells
Two daughter cells will be in the same place of the adult cell
"""
function divide(phyt)
    phytos = zeros(Real,size(phyt,1)*2)
    # NOT all C and N can turn into new cells
    for i in 1:2
        phytos[1+(i-1)*12]  = phyt[1]          # x
        phytos[2+(i-1)*12]  = phyt[2]          # y
        phytos[3+(i-1)*12]  = phyt[3]          # z
        phytos[4+(i-1)*12]  = phyt[4]          # species
        phytos[5+(i-1)*12]  = phyt[5] .+ 1.0   # generation
        phytos[6+(i-1)*12]  = 1.0              # age
        phytos[7+(i-1)*12]  = phyt[7] .* 0.45  # size
        phytos[8+(i-1)*12]  = phyt[8] .* 0.5   # Cq1
        phytos[9+(i-1)*12]  = phyt[9] .* 0.45  # Cq2
        phytos[10+(i-1)*12] = phyt[10] .* 0.5  # Nq
        phytos[11+(i-1)*12] = phyt[11] .* 0.5  # Pq
        phytos[12+(i-1)*12] = phyt[12] .* 0.5  # chl
    end
    return phytos
end
###################################
# model update for phytoplanktons #
###################################
"""
    phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
Update the individuals of current time step into next time step
Control flow: Graze ->Grow(photosynthesis, biosynthesis, maintenance) -> Natural Death -> Division
Return a dataframe of next time step individuals, graze number, divide number, death number, and nutrient consumption
"""
function phyt_update(model, ΔT::Int64)
    t = model.t
    g = model.grid
    nutrients = model.nutrients
    IR = model.PAR
    temp = model.temp
    params = model.params
    phyts_a = copy(model.individuals.phytos)

    # load nutrients
    counts = pop_counts()
    chl_num = count_chl(phyts_a, g)
    cumsum_chl = cumsum(chl_num, dims = 3)

    #set up a empty array to record all updated agents
    phyts_b = Real[]
    consume = nutrients_init(g)
    # iterate phytoplankton agents
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        sp = Int(phyt[4])
        x, y, z = which_grid(phyt, g)
        temp_t = temp[x,y,z,t]
        IR_t = IR[x,y,t]
        NH4 = max(0.0, nutrients.NH4[x, y, z])
        NO3 = max(0.0, nutrients.NO3[x, y, z])
        PO4 = max(0.0, nutrients.PO4[x, y, z])

        # compute probabilities of grazing
        # Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        if params["Grz_P"] == 0
            P_graz = false
        else
            reg_graz = phyt[7]/params["Grz_P"]
            P_graz = rand(Bernoulli(reg_graz))
        end

        if P_graz == false #not grazed
            # compute death probability after a certain age, may be abandoned
            reg_age = max(0.0, phyt[6] - params["death_age"])
            shape_factor_death = params["a_death"]*reg_age^params["b_death"]
            P_death = rand(Bernoulli(shape_factor_death/(1+shape_factor_death)))

            if P_death == false # not natural death
                # Compute light attenuation
                atten = (params["katten_w"] + params["katten_c"] * cumsum_chl[x, y, z])*(-g.zF[z])
                if atten == 0.0
                    α_I = params["α"]*IR_t
                else
                    α_I = params["α"]*IR_t*(1.0 - exp(-atten))/atten
                end

                # Compute photosynthesis rate
                Tempstd = exp(params["TempAe"]*(1.0/(temp_t+273.15)-1.0/params["Tempref"]))
                photoTempFunc = params["TempCoeff"]*max(1.0e-10,Tempstd)
                PCmax_sp = params["PCmax"][sp]/86400
                PCm = PCmax_sp*photoTempFunc*phyt[7]^params["PC_b"][sp]
                PC = PCm*(1-exp(-α_I*phyt[12]/(phyt[8]*PCm)))
                PS = PC*phyt[8] # unit: mmol C/second/individual
                Eₖ = PCm/(phyt[12]/phyt[8]*params["α"])
                tmp = α_I/params["α"]
                if (tmp > Eₖ) & (params["inhibcoef"][sp] == 1.0)
                    PS = PS*Eₖ/tmp*params["inhibcoef"][sp]
                end
                PP = PS*ΔT # unit: mmol C/hour/individual

                # Compute cell-based N uptake rate according Droop limitation
                Qn = (phyt[10]+phyt[8]*params["R_NC"])/(phyt[8]+phyt[9])
                #In-Cell N uptake limitation
                regQn = max(0.0,min(1.0,(params["Nqmax"]-Qn)/(params["Nqmax"]-params["Nqmin"])))
                VNH4max_sp = params["VNH4max"][sp]/86400
                VNO3max_sp = params["VNO3max"][sp]/86400
                VNH4m = VNH4max_sp*phyt[7]^params["VN_b"][sp]
                VNO3m = VNO3max_sp*phyt[7]^params["VN_b"][sp]
                NH4uptake = VNH4m*NH4/(NH4+params["KsatNH4"][sp])*regQn
                NO3uptake = VNO3m*NO3/(NO3+params["KsatNO3"][sp])*regQn
                VNH4cell = NH4uptake*phyt[8] # unit: mmol N/second/individual
                VNO3cell = NO3uptake*phyt[8] # unit: mmol N/second/individual
                VNH4 = min(NH4*g.V[x,y,z]/10.0, VNH4cell*ΔT) # unit: mmol N/hour/individual
                VNO3 = min(NO3*g.V[x,y,z]/10.0, VNO3cell*ΔT) # unit: mmol N/hour/individual

                # Compute cell-based P uptake rate according Droop limitation
                Qp = (phyt[11]+phyt[8]*params["R_PC"])/(phyt[8]+phyt[9])
                #In-Cell P uptake limitation
                regQp = max(0.0,min(1.0,(params["Pqmax"]-Qp)/(params["Pqmax"]-params["Pqmin"])))
                VPmax_sp = params["VPmax"][sp]/86400
                VPm = VPmax_sp*phyt[7]^params["VP_b"][sp]
                Puptake = VPm*PO4/(PO4+params["KsatP"][sp])*regQp
                VPcell = Puptake*phyt[8] # unit: mmol P/second/individual
                VPO4 = min(PO4*g.V[x,y,z]/10.0, VPcell*ΔT) # unit: mmol P/hour/individual

                # Compute the ratio of chl synthesis and N uptake
                # ρ equals to ratio of the realised quantum efficiency for photosynthesis divided by the maximum efficiency
                if IR_t > 0
                    ρ_chl = PS*params["Chl2N"]/(params["α"]*IR_t*phyt[12]/phyt[8])
                else
                    ρ_chl = 0.0
                end

                # Metabolic partitioning for biosynthesis, decrease with size
                shape_factor_β = params["a_β"]*phyt[7]^params["b_β"]
                β = shape_factor_β/(1+shape_factor_β)

                # Compute extra cost for biosynthesis, return a rate (per hour)
                respir_extra = params["respir_ex"]*phyt[7]^params["respir_b"]

                # C, N, P storages update
                phyt[9] = phyt[9] + PP
                phyt[10]= phyt[10]+ VNH4 + VNO3
                phyt[11]= phyt[11]+ VPO4

                # maximum biosynthesis rate
                SynC = β*params["k_mtb"]*phyt[9]/(1+respir_extra)

                # maximum allowed biosynthesis rate by Nq and Pq
                BrNC = phyt[10]/params["R_NC"]
                BrPC = phyt[11]/params["R_PC"]
                Bmax = min(BrNC, BrPC)

                # excretion
                if SynC > Bmax
                    excretC = SynC - Bmax
                    SynC = Bmax
                else
                    excretC = 0.0
                end

                # update quotas and biomass
                CostC = SynC*(1+respir_extra)
                MaintenC = (1-β)*SynC/β
                phyt[8] = phyt[8] + SynC
                phyt[9] = phyt[9] - CostC -MaintenC
                phyt[10]= phyt[10]- SynC*params["R_NC"]
                phyt[11]= phyt[11]- SynC*params["R_PC"]

                dsize= SynC/(phyt[8])
                phyt[7]  = max(0.0,phyt[7]+dsize)
                phyt[12] = phyt[12] + ρ_chl*SynC*params["R_NC"]
                phyt[6]  = phyt[6] + 1.0*(ΔT/3600)

                # compute probabilities of division
                reg_size = params["dvid_stp"]*(phyt[7] - params["dvid_size"])
                reg_divide = 0.05*(tanh(reg_size) + 1)
                P_dvi = rand(Bernoulli(reg_divide))
                if P_dvi == false # not divide
                    append!(phyts_b,phyt)
                else # divide
                    counts.divid += 1
                    phyts = divide(phyt)
                    append!(phyts_b,phyts)
                    consume.DIC[x, y, z] = consume.DIC[x, y, z] + phyt[9]*0.1 # consume C when cell is divided
                end # divide
                consume.DIC[x, y, z] = consume.DIC[x, y, z] + MaintenC + CostC - SynC
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + excretC
                consume.NH4[x, y, z] = consume.NH4[x, y, z] - VNH4
                consume.NO3[x, y, z] = consume.NO3[x, y, z] - VNO3
                consume.PO4[x, y, z] = consume.PO4[x, y, z] - VPO4
            else # natural death
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt[8]+phyt[9])*params["mortFracC"]
                consume.DON[x, y, z] = consume.DON[x, y, z] + phyt[10]*params["mortFracN"]
                consume.DOP[x, y, z] = consume.DOP[x, y, z] + phyt[11]*params["mortFracP"]
                consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt[8]+phyt[9])*(1.0 - params["mortFracC"])
                consume.PON[x, y, z] = consume.PON[x, y, z] + phyt[10]*(1.0 - params["mortFracN"])
                consume.POP[x, y, z] = consume.POP[x, y, z] + phyt[11]*(1.0 - params["mortFracP"])
                counts.death += 1
            end # naturan death
        else #grazed, no sloppy feeding here, all nutrients go back to organic pools
            counts.graze += 1
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt[8]+phyt[9])*params["grazFracC"]
            consume.DON[x, y, z] = consume.DON[x, y, z] + phyt[10]*params["grazFracN"]
            consume.DOP[x, y, z] = consume.DOP[x, y, z] + phyt[11]*params["grazFracP"]
            consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt[8]+phyt[9])*(1.0 - params["grazFracC"])
            consume.PON[x, y, z] = consume.PON[x, y, z] + phyt[10]*(1.0 - params["grazFracN"])
            consume.POP[x, y, z] = consume.POP[x, y, z] + phyt[11]*(1.0 - params["grazFracP"])
        end # graze
    end # while loop to traverse the array of agents
    phyts_b = reshape(phyts_b,size(phyts_a,1),Int(length(phyts_b)/size(phyts_a,1)))
    return phyts_b,counts,consume
end # for loop of time
