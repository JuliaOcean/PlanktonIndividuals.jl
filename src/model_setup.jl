"""
    setup_agents(RunParams,grid)
Set up a series of agents following a normal distribution (mean,var)
'Nindivi' is agent number for each species, 'sp' is number of species, 'Nsuper' is the number of cells one agent represents,
'Cquota' is the initial biomass for one cell
'(x,y,z)' of an individual is the actual location not grid indices
"""
function setup_agents(RunParam::RunParams,grid)
    PhytoOpt = RunParam.PhytoOpt
    phyts0 = DataFrame(x=Float64[], y=Float64[], z=Float64[],
                       gen=Int64[], size=Float64[], Cq1=Float64[],
                       Cq2=Float64[], Nq=Float64[], Pq=Float64[],
                       chl=Float64[], sp=Int64[], age=[])
    for i in 1:PhytoOpt.Nsp
        for j in 1:PhytoOpt.Nindivi
            # agent location (actual location)
            x = rand(Uniform(grid.xF[1],grid.xF[end]))
            y = rand(Uniform(grid.yF[1],grid.yF[end]))
            z = rand(Uniform(grid.zF[1],grid.zF[end]))
            # a normal distribution with mean variance
            radm = max(0.05, rand(Normal(PhytoOpt.mean,PhytoOpt.var)))
            gen  = 1
            size = radm
            Cq1  = PhytoOpt.Cquota[i]*PhytoOpt.Nsuper # Nsuper is the number of cells one super agent repersents
            Cq2  = PhytoOpt.Cquota[i]*PhytoOpt.Nsuper*radm
            Nq   = 13/120*Cq2
            Pq   = 1/120*Cq2
            chl  = Cq2*0.4 # mgChl(/mmolC)
            sp   = i
            age  = 1.0
            push!(phyts0,(x=x,y=y,z=z,gen=gen,size=size,
                          Cq1=Cq1,Cq2=Cq2,Nq=Nq,Pq=Pq,
                          chl=chl,sp=sp,age=age))
        end
    end
    B = [phyts0]
    if RunParam.Zoo == false
        return B
    else
        ZP = setup_zooplkt(RunParam.ZooOpt, grid)
        B = hcat(B,ZP)
        return B
    end
end

"""
    setup_zooplkt(ZooOpt, grid)
Set up zooplankton individuals according to 'ZooOpt' from 'RunParam'
"""
function setup_zooplkt(ZooOpt, grid)
    zoo0 = DataFrame(x=Float64[], y=Float64[], z=Float64[],
                     gen=Int64[], size=Float64[], Cq1=Float64[],
                     Cq2=Float64[], Nq=Float64[], Pq=Float64[],
                     chl=Float64[], sp=Int64[], age=[])
    for i in 1:ZooOpt.Nsp
        for j in 1:ZooOpt.Nindivi
            # agent location (actual location)
            x = rand(Uniform(grid.xF[1],grid.xF[end]))
            y = rand(Uniform(grid.yF[1],grid.yF[end]))
            z = rand(Uniform(grid.zF[1],grid.zF[end]))
            # a normal distribution with mean variance
            radm = max(0.05, rand(Normal(ZooOpt.mean,ZooOpt.var)))
            gen  = 1
            size = radm
            Cq1  = 0.0  # zooplankton has only one C quota
            Cq2  = ZooOpt.Cquota[i]*ZooOpt.Nsuper*radm
            Nq   = 13/120*Cq2
            Pq   = 1/120*Cq2
            chl  = 0.0 # no Chl in zooplankton
            sp   = i
            age  = 1.0
            push!(zoo0,(x=x,y=y,z=z,gen=gen,size=size,
                        Cq1=Cq1,Cq2=Cq2,Nq=Nq,Pq=Pq,
                        chl=chl,sp=sp,age=age))
        end
    end
    return [zoo0]
end


"""
    nutrients_init(g)
"""
function nutrients_init(g)
    nut = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    return nut
end


"""
    setup_nutrients(g,nut)
Set up initial nutrient fields according to grid information
'Nut' is an array of 10 elements, each element is a kind of nutrient
"""
function setup_nutrients(g,nut)
    DIC = fill(nut[1],(g.Nx, g.Ny, g.Nz))
    NH4 = fill(nut[2],(g.Nx, g.Ny, g.Nz))
    NO3 = fill(nut[3],(g.Nx, g.Ny, g.Nz))
    PO4 = fill(nut[4],(g.Nx, g.Ny, g.Nz))
    DOC = fill(nut[5],(g.Nx, g.Ny, g.Nz))
    DON = fill(nut[6],(g.Nx, g.Ny, g.Nz))
    DOP = fill(nut[7],(g.Nx, g.Ny, g.Nz))
    POC = fill(nut[8],(g.Nx, g.Ny, g.Nz))
    PON = fill(nut[9],(g.Nx, g.Ny, g.Nz))
    POP = fill(nut[10],(g.Nx, g.Ny, g.Nz))

    nutrients = nutrient_fields(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP)
    return nutrients
end


