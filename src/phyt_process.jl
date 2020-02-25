###########################################
# define functions of pysiology processes #
###########################################
"""
    divide(phyt)
An adult cell divides evenly into two daughter cells
Two daughter cells will be in the same place of the adult cell
"""
function divide(phyt::DataFrameRow)
    phytops = DataFrame(x=Float64[0.0,0.0], y=Float64[0.0,0.0], z=Float64[0.0,0.0], gen=Int64[1,1], size=Float64[0.0,0.0], Cq1=Float64[0.0,0.0], Cq2=Float64[0.0,0.0], Nq=Float64[0.0,0.0], chl=Float64[0.0,0.0], sp=Int64[0,0], age=Float64[0.0,0.0])  # initialize new cell
    # NOT all C and N can turn into new cells
    phytops[1,:].x = phyt.x
    phytops[1,:].y = phyt.y
    phytops[1,:].z = phyt.z
    phytops[1,:].gen = phyt.gen + 1
    phytops[1,:].Cq1 = phyt.Cq1 * 0.5
    phytops[1,:].Cq2 = phyt.Cq2 * 0.45
    phytops[1,:].Nq  = phyt.Nq  * 0.5
    phytops[1,:].size= phyt.size* 0.45
    phytops[1,:].chl = phyt.chl * 0.5
    phytops[1,:].sp = phyt.sp
    phytops[1,:].age = 1.0

    phytops[2,:].x = phyt.x
    phytops[2,:].y = phyt.y
    phytops[2,:].z = phyt.z
    phytops[2,:].gen = phyt.gen + 1
    phytops[2,:].Cq1 = phyt.Cq1 * 0.5
    phytops[2,:].Cq2 = phyt.Cq2 * 0.45
    phytops[2,:].Nq  = phyt.Nq  * 0.5
    phytops[2,:].size= phyt.size* 0.45
    phytops[2,:].chl = phyt.chl * 0.5
    phytops[2,:].sp = phyt.sp
    phytops[2,:].age = 1.0

    return phytops
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
function phyt_update(t::Int64, ΔT::Int64, phyts_a, model)
    g = model.grid
    nutrients = model.nutrients
    IR = model.PAR
    temp = model.temp
    params = model.params

    # load nutrients
    dvid_ct = 0; graz_ct = 0; death_ct = 0
    Num_phyt = size(phyts_a,1)
    chl_num = count_chl(phyts_a, g)
    cumsum_chl = cumsum(chl_num, dims = 3)

    #set up a dataframe to record all updated agents
    phyts_b = DataFrame(x=Float64[], y=Float64[], z=Float64[],
                        gen=Int64[], size=Float64[], Cq1=Float64[],
                        Cq2=Float64[], Nq=Float64[], chl=Float64[],
                        sp=Int64[], age=Float64[])

    consume = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                              zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                              zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    # iterate phytoplankton agents
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        sp = phyt.sp
        z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
        temp_t = temp[x,y,z,t]
        IR_t = IR[x,y,t]
        DIN = max(0.0, nutrients.DIN[x, y, z])

        # compute probabilities of grazing
        # Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        if params["Grz_P"] == 0
            P_graz = false
        else
            reg_num_graz = exp(Num_phyt/RunParam.PhytoOpt.Nindivi/RunParam.PhytoOpt.Nsp)
            P_graz = rand(Bernoulli(reg_num_graz*phyt.size/params["Grz_P"]))
        end

        # compute death probability after a certain age, may be abandoned
        reg_age = max(0.0, phyt.age - params["death_age"])
        shape_factor_death = params["a_death"]*reg_age^params["b_death"]
        P_death = rand(Bernoulli(shape_factor_death/(1+shape_factor_death)))

        # Compute extra cost for biosynthesis, return a rate (per hour)
        respir_extra = params["respir_ex"]*phyt.size^params["respir_b"]

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
        PCmax_sp = params["PCmax"][phyt.sp]/86400
        PCm = PCmax_sp*photoTempFunc*phyt.size^params["PC_b"][phyt.sp]
        PC = PCm*(1-exp(-α_I*phyt.chl/(phyt.Cq2*PCm)))
        PS = PC*phyt.Cq2 # unit: mmol C/second/individual
        Eₖ = PCm/(phyt.chl/phyt.Cq2*params["α"])
        tmp = α_I/params["α"]
        if (tmp > Eₖ) & (params["inhibcoef"][phyt.sp] == 1.0)
            PS = PS*Eₖ/tmp*params["inhibcoef"][phyt.sp]
        end
        PP = PS*ΔT # unit: mmol C/hour/individual

        # Compute cell-based N uptake rate according Droop limitation
        Qn = phyt.Nq/(phyt.Cq1+phyt.Cq2)
        #In-Cell N uptake limitation
        regQ = max(0.0,min(1.0,(params["Nqmax"]-Qn)/(params["Nqmax"]-params["Nqmin"])))
        VNmax_sp = params["VNmax"][phyt.sp]/86400
        VNm = VNmax_sp*phyt.size^params["VN_b"][phyt.sp]
        Nuptake = VNm*DIN/(DIN+params["KsatN"][phyt.sp])*regQ
        VNcell = Nuptake*phyt.Cq2 # unit: mmol N/second/individual
        VN = min(DIN*g.V[x,y,z]/10.0, VNcell*ΔT) # unit: mmol N/hour/individual

        # Compute the ratio of chl synthesis and N uptake
        # ρ equals to ratio of the realised quantum efficiency for photosynthesis divided by the maximum efficiency
        if IR_t > 0
            ρ_chl = PS/(params["α"]*IR_t*phyt.chl/phyt.Cq2)
        else
            ρ_chl = 0.0
        end

        # Metabolic partitioning for biosynthesis, decrease with size
        shape_factor_β = params["a_β"]*phyt.size^params["b_β"]
        β = shape_factor_β/(1+shape_factor_β)

        if P_graz == false #not grazed
            SynC = min(VN/params["R_NC"],β*params["k_mtb"]*phyt.Cq1/(1+respir_extra))
            CostC = SynC*(1+respir_extra)
            MaintenC = (1-β)*params["k_mtb"]*phyt.Cq1
            dCq1 = PP - CostC - MaintenC
            dCq2 = SynC
            dNq  = VN
            dsize= dCq2/(phyt.Cq2)
            phyt.Cq1 = max(0.0, phyt.Cq1 + dCq1)
            phyt.Cq2 = phyt.Cq2 + dCq2
            phyt.Nq  = phyt.Nq + dNq
            phyt.size= max(0.0,phyt.size+dsize)
            phyt.chl = phyt.chl + ρ_chl*VN*params["Chl2N"]
            phyt.age = phyt.age + 1.0*(ΔT/3600)
            if P_death == false # not natural death
                # compute probabilities of division
                reg_size = max(0.0, phyt.size - params["dvid_size"])
                shape_factor_divide = (params["a_dvi"][sp]*reg_size)^params["b_dvi"][sp]
                P_dvi = rand(Bernoulli(shape_factor_divide/(1+shape_factor_divide)))
                if P_dvi == false # not divide
                    push!(phyts_b,phyt)
                else # divide
                    dvid_ct += 2
                    phyts2 = divide(phyt)
                    append!(phyts_b,phyts2)
                    consume.DIC[x, y, z] = consume.DIC[x, y, z] + phyt.Cq2*0.1 # consume C when cell is divided
                end # divide
            else # natural death
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt.Cq1+phyt.Cq2)*params["mortFracC"]
                consume.DON[x, y, z] = consume.DON[x, y, z] + phyt.Nq*params["mortFracN"]
                consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt.Cq1+phyt.Cq2)*(1.0 - params["mortFracC"])
                consume.PON[x, y, z] = consume.PON[x, y, z] + phyt.Nq*(1.0 - params["mortFracN"])
                death_ct += 1
            end # naturan death
            consume.DIC[x, y, z] = consume.DIC[x, y, z] + MaintenC + CostC - SynC
            consume.DIN[x, y, z] = consume.DIN[x, y, z] - VN
        else #grazed
            graz_ct += 1
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt.Cq1+phyt.Cq2)*params["grazFracC"]*0.5
            consume.DON[x, y, z] = consume.DON[x, y, z] + phyt.Nq*params["grazFracN"]*0.5
            consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt.Cq1+phyt.Cq2)*(1.0 - params["grazFracC"])*0.5
            consume.PON[x, y, z] = consume.PON[x, y, z] + phyt.Nq*(1.0 - params["grazFracN"])*0.5
        end # graze
    end # while loop to traverse the array of agents
    return phyts_b,(dvid_ct,graz_ct,death_ct),consume
end # for loop of time
