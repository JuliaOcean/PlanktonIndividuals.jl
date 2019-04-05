###########################################
# define functions of pysiology processes #
###########################################
k_respir(Cq2) = respir_a*(12.0e9*Cq2)^respir_b/Cq2
respir(Cq2) = respir_a*(12.0e9*Cq2)^respir_b

function daynight(t::Int64, IR)
    if IR[trunc(Int,t*deltaT/3600)] < 5.0
        return false
    else
        return true
    end
end

function PAR_cal(I, z, cumsum_cell)
    atten = (katten_w *(-z) + katten_c * cumsum_cell)
    PAR = α*I*exp(-atten)
    return PAR
end

function PC(PAR, Temp, phyt) 
    Tempstd = exp(TempAe*(1.0/(Temp+273.15)-1.0/Tempref))
    photoTempFunc = TempCoeff*max(1.0e-10,Tempstd)
    PC = PCmax[phyt.sp]*photoTempFunc*(1-exp(-PAR*phyt.chl/phyt.Cq2/PCmax[phyt.sp]))*Cquota[phyt.sp]*phyt.size
    return PC
end

function Nuptake(Nit, phyt)
    Nqmax = Nqmax_a*phyt.size*Cquota[phyt.sp]
    Nqmin = Nqmin_a*phyt.size*Cquota[phyt.sp]
    #In-Cell N uptake limitation
    regQ = max(0.0,min(1.0,(Nqmax-phyt.Nq)/(Nqmax-Nqmin)))
    Nuptake = VNmax[phyt.sp]*Nit/(Nit+KsatN)*regQ*Cquota[phyt.sp]*phyt.size
    return Nuptake
end

function chl_sync(phyt,PP,I)
    if I > 0
        ρ_chl = PP/(α*I*phyt.chl)
    else
        ρ_chl = 1.0
    end
    return ρ_chl
end

function divide(phyt::DataFrameRow)
    phytops = DataFrame(x=Float64[0.0,0.0], y=Float64[0.0,0.0], z=Float64[0.0,0.0], gen=Int64[1,1], size=Float64[0.0,0.0], Cq1=Float64[0.0,0.0], Cq2=Float64[0.0,0.0], Nq=Float64[0.0,0.0], chl=Float64[0.0,0.0], sp=Int64[0,0])  # initialize new cell
    # NOT all C and N can turn into new cells
    phytops[1,:].x = phyt.x
    phytops[1,:].y = phyt.y
    phytops[1,:].z = phyt.z
    phytops[1,:].gen = phyt.gen + 1
    phytops[1,:].Cq1 = phyt.Cq1 * 0.5
    phytops[1,:].Cq2 = phyt.Cq2 * 0.5
    phytops[1,:].Nq  = phyt.Nq  * 0.5
    phytops[1,:].size= phyt.size* 0.5
    phytops[1,:].chl = phyt.chl * 0.5
    phytops[1,:].sp = phyt.sp

    phytops[2,:].x = phyt.x
    phytops[2,:].y = phyt.y
    phytops[2,:].z = phyt.z
    phytops[2,:].gen = phyt.gen + 1
    phytops[2,:].Cq1 = phyt.Cq1 * 0.5
    phytops[2,:].Cq2 = phyt.Cq2 * 0.5
    phytops[2,:].Nq  = phyt.Nq  * 0.5
    phytops[2,:].size= phyt.size* 0.5
    phytops[2,:].chl = phyt.chl * 0.5
    phytops[2,:].sp = phyt.sp

    return phytops
end
####################################
# model update (ONE time step: 1h) #
####################################
function update(t::Int64, deltaT::Int64, phyts_a, nutrients, IR, temp, cell_num)
    # load nutrients
    DIN = copy(nutrients.DIN[t])
    DON = copy(nutrients.DON[t])
    PON = copy(nutrients.PON[t])
    DOC = copy(nutrients.DOC[t])
    POC = copy(nutrients.POC[t])
    dvid_ct  = 0
    graz_ct  = 0
    death_ct = 0
    Num_phyt = size(phyts_a,1)
    cumsum_cell = cumsum(cell_num, dims = 3)
    #set up a dataframe to record all updated agents
    phyts_b = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[])
    #
    # iterate phytoplankton agents
    #
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
        #compute probabilities of grazing and division
        P_graz = rand(Bernoulli(exp(Num_phyt/N*Nsp)*phyt.size/Grz_P))
        # Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        P_dvi=max(0.0,phyt.size-dvid_size)*1.0e5*rand(Bernoulli(phyt.size/Dvid_P))
        PAR = PAR_cal(IR[trunc(Int,t*deltaT/3600)], g.zF[z], cumsum_cell[y, x, z])
        PP = PC(PAR,temp[trunc(Int,t*deltaT/3600)],phyt)*deltaT
        VN = Nuptake(DIN,phyt)*deltaT
        Dmd_NC = (1+k_respir(phyt.Cq1))*VN/R_NC
        Res2 = respir(phyt.Cq2)*deltaT
        ρ_chl = chl_sync(phyt,PP,IR[trunc(Int,t*deltaT/3600)])
        if P_graz < 1 #not grazed
            if phyt.Cq2+phyt.Cq1 ≤ Cmin # natural death
                DOC = DOC + (phyt.Cq1+phyt.Cq2)*mortFracC
                DON = DON + phyt.Nq*mortFracN
                POC = POC + (phyt.Cq1+phyt.Cq2)*(1.0 - mortFracC)
                PON = PON + phyt.Nq*(1.0 - mortFracN)
                death_ct += 1
            else # not natural death
                if daynight(t,IR) # day
                    if PP > Dmd_NC
                        ExuC = (PP - Dmd_NC)*FracExuC
                        CostC= Dmd_NC
                        SynC = VN/R_NC
                    else
                        CostC= PP
                        SynC = PP/(1+k_respir(phyt.Cq1))
                        ExuC = 0.0
                        VN   = SynC*R_NC
                    end #exudation
                    dCq1 = PP - CostC - ExuC
                    dCq2 = SynC - Res2
                    dNq  = VN - Res2*R_NC
                    dsize= max((dCq1+dCq2)/(phyt.Cq1+phyt.Cq2),dNq/phyt.Nq)
                    phyt.Cq1 = max(0.0,phyt.Cq1 + dCq1)
                    phyt.Cq2 = max(0.0,phyt.Cq2 + dCq2)
                    phyt.Nq  = max(0.0,phyt.Nq + dNq)
                    phyt.size= phyt.size+dsize
                    phyt.chl = phyt.chl + ρ_chl*VN*Chl2N
                    push!(phyts_b,phyt)
                    DIN = DIN - VN
                    DOC = DOC + ExuC
                else # night
                    if P_dvi < 1 # not divide
                        CostC= min(0.4*phyt.Cq1,Dmd_NC)
                        SynC = min(VN/R_NC,0.4*phyt.Cq1/(1+k_respir(phyt.Cq1)))
                        ExuC = 0.0
                        VN   = SynC*R_NC
                        dCq1 = PP - CostC - ExuC
                        dCq2 = SynC - Res2
                        dNq  = VN - Res2*R_NC
                        dsize= max((dCq1+dCq2)/(phyt.Cq1+phyt.Cq2),dNq/phyt.Nq)
                        phyt.Cq1 = max(0.0,phyt.Cq1 + dCq1)
                        phyt.Cq2 = max(0.0,phyt.Cq2 + dCq2)
                        phyt.Nq  = max(0.0,phyt.Nq + dNq)
                        phyt.size= phyt.size+dsize
                        phyt.chl = phyt.chl + ρ_chl*VN*Chl2N
                        push!(phyts_b,phyt)
                        DIN = DIN - VN
                        DOC = DOC + ExuC
                    else #divide
                        dvid_ct += 2
                        global dvdcount += 2
                        phyts2 = divide(phyt)
                        append!(phyts_b,phyts2)
                        DOC = DOC + phyt.Cq2*0.0 # consume C when cell is divided
                    end # divide
                end # day night?
            end # natural death
        else #grazed
            graz_ct += 1
            DOC = DOC + (phyt.Cq1+phyt.Cq2)*grazFracC*0.5
            DON = DON + phyt.Nq*grazFracN*0.5
            POC = POC + (phyt.Cq1+phyt.Cq2)*(1.0 - grazFracC)*0.5
            PON = PON + phyt.Nq*(1.0 - grazFracN)*0.5
        end # graze
    end # while loop to traverse the array of agents
    #
    # update nutrients
    sinkPOC  = POC*k_sink*deltaT
    sinkPON  = PON*k_sink*deltaT
    reminPOC = POC*kPOC*deltaT
    reminPON = PON*kPON*deltaT
    reminDOC = DOC*kDOC*deltaT
    reminDON = DON*kDON*deltaT

    PON = PON - sinkPON - reminPON
    POC = POC - sinkPOC - reminPOC
    DON = DON + reminPON - reminDON
    DOC = DOC + reminPOC - reminDOC
    DIN = DIN + reminDON
    push!(nutrients,(DIN=DIN, DOC=DOC, DON=DON, POC=POC, PON=PON))
    return phyts_b,dvid_ct,graz_ct
end # for loop of time
