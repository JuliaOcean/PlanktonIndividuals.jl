###########################################
# define functions of pysiology processes #
###########################################
k_respir(Cq2) = respir_a*(12.0e9*Cq2)^respir_b/Cq2
respir(Cq2) = respir_a*(12.0e9*Cq2)^respir_b

function daynight(t::Int64, IR)
    if IR[t] < 5.0
        return false
    else
        return true
    end
end

function PAR_cal(phyt, I, zf, cell_num)
    cumsum_cell = cumsum(cell_num, dims = 3)
    z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
    atten = (katten_w *(-zf[z]) + katten_c * cumsum_cell[y, x, z])
    PAR = α*I*exp(-atten)
    return PAR
end

function PC(PAR, Temp, phyt) 
    Tempstd = exp(TempAe*(1.0/(Temp+273.15)-1.0/Tempref))
    photoTempFunc = TempCoeff*max(1.0e-10,Tempstd)
    PC = PCmax*photoTempFunc*(1-exp(-PAR*phyt.chl/phyt.Cq2/PCmax))*Cquota[phyt.sp]*phyt.size
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
    phytops[1,:].Cq1 = phyt.Cq1 * 0.45
    phytops[1,:].Cq2 = phyt.Cq2 * 0.45
    phytops[1,:].Nq  = phyt.Nq  * 0.45
    phytops[1,:].size= phyt.size* 0.45
    phytops[1,:].chl = phyt.chl * 0.45
    phytops[1,:].sp = phyt.sp

    phytops[2,:].x = phyt.x
    phytops[2,:].y = phyt.y
    phytops[2,:].z = phyt.z
    phytops[2,:].gen = phyt.gen + 1
    phytops[2,:].Cq1 = phyt.Cq1 * 0.45
    phytops[2,:].Cq2 = phyt.Cq2 * 0.45
    phytops[2,:].Nq  = phyt.Nq  * 0.45
    phytops[2,:].size= phyt.size* 0.45
    phytops[2,:].chl = phyt.chl * 0.45
    phytops[2,:].sp = phyt.sp

    return phytops
end
######################################
# advection and convection of agents #
######################################
function trilinear_itpl(x,y,z,vel_field,t::Int64)
    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    xᵈ = x - x₀
    yᵈ = y - y₀
    zᵈ = z - z₀

    vel_000 = vel_field[y₀, x₀, z₀, t]
    vel_100 = vel_field[y₀+1, x₀, z₀, t]
    vel_001 = vel_field[y₀, x₀, z₀+1, t]
    vel_010 = vel_field[y₀, x₀+1, z₀, t]
    vel_110 = vel_field[y₀+1, x₀+1, z₀, t]
    vel_011 = vel_field[y₀, x₀+1, z₀+1, t]
    vel_101 = vel_field[y₀+1, x₀, z₀+1, t]
    vel_111 = vel_field[y₀+1, x₀+1, z₀+1, t]

    vel_00 = vel_000 * (1 - yᵈ) + vel_100 * yᵈ
    vel_01 = vel_001 * (1 - yᵈ) + vel_101 * yᵈ
    vel_10 = vel_010 * (1 - yᵈ) + vel_110 * yᵈ
    vel_11 = vel_011 * (1 - yᵈ) + vel_111 * yᵈ

    vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
    vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ

    vel = vel_0 * (1-zᵈ) + vel_1 * zᵈ

    return vel
end

function agent_move(phyts_a,bdry,u,v,w,xgrid,ygrid,zgrid,t::Int64)
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        uvel = trilinear_itpl(phyt.x, phyt.y, phyt.z, u, t)
        vvel = trilinear_itpl(phyt.x, phyt.y, phyt.z, v, t)
        wvel = trilinear_itpl(phyt.x, phyt.y, phyt.z, w, t)

        xi, yi, zi = trunc(Int,phyt.x), trunc(Int,phyt.y), trunc(Int,phyt.z)
        dx = uvel*3600/1000/xgrid[xi]/96.4 # 1 degree of lat at 30N
        dy = vvel*3600/1000/ygrid[yi]/111 # 1 degree of lon
        dz = wvel*3600/zgrid[zi]
        phyt.x = phyt.x - dx*(1+rand()/5)
        phyt.y = phyt.y - dy*(1+rand()/5)
        phyt.z = max(bdry[3,1],min(bdry[3,2],phyt.z - dz*(1+rand()/5)))
	# periodic domian
	if phyt.x ≥ bdry[1,2]
		phyt.x = phyt.x - bdry[1,2]
	end
	if phyt.x ≤ bdry[1,1]
		phyt.x = phyt.x + bdry[1,2]
	end
	if phyt.y ≥ bdry[2,2]
		phyt.y = phyt.y - bdry[2,2]
	end
	if phyt.y ≤ bdry[2,1]
		phyt.y = phyt.y + bdry[2,2]
	end
    end
end
####################################
# model update (ONE time step: 1h) #
####################################
function update(t::Int64, phyts_a, nutrients, IR, temp, cell_num)
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
    #set up a dataframe to record all updated agents
    phyts_b = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[])
#
# iterate phytoplankton agents
#
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        #compute probabilities of grazing and division
	P_graz = rand(Bernoulli(exp(Num_phyt/N*Nsp)*phyt.size/Grz_P))
# Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        P_dvi=max(0.0,phyt.size-dvid_size)*1.0e5*rand(Bernoulli(phyt.size/Dvid_P))
	PAR = PAR_cal(phyt, IR[t], zf, cell_num)
        PP = PC(PAR,temp[t],phyt)
        VN = Nuptake(DIN,phyt)
        Dmd_NC = (1+k_respir(phyt.Cq1))*VN/R_NC
        Res2 = respir(phyt.Cq2)
        ρ_chl = chl_sync(phyt,PP,IR[t])
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
    sinkPOC  = POC*k_sink
    sinkPON  = PON*k_sink
    reminPOC = POC*kPOC
    reminPON = PON*kPON
    reminDOC = DOC*kDOC
    reminDON = DON*kDON

    PON = PON - sinkPON - reminPON
    POC = POC - sinkPOC - reminPOC
    DON = DON + reminPON - reminDON
    DOC = DOC + reminPOC - reminDOC
    DIN = DIN + reminDON
    push!(nutrients,(DIN=DIN, DOC=DOC, DON=DON, POC=POC, PON=PON))
    return phyts_b,dvid_ct,graz_ct
end # for loop of time