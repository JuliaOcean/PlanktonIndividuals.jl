# set up a series of agents following a normal distribution
function setup_agents(N::Int64,Cquota::Array,mean::Float64,var::Float64,grid)
    phyts0 = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[])
    for i in 1:N
        # agent location
        x = rand(1.5*10:(grid.Nx-0.5)*10)/10
        y = rand(1.5*10:grid.Ny*5)/10
        z = rand(1.5*10:grid.Nz*7.5)/10
        # a normal distribution with mean variance
        radm = max(0.05, rand(Normal(mean,var)))
        gen  = 1
        size = radm
        Cq1  = Cquota[1]*1E-3
        Cq2  = Cquota[1]*radm
        Nq   = 13/106*2*Cq2
        chl  = Cq2*0.4
        sp   = 1
        push!(phyts0,(x=x,y=y,z=z,gen=gen,size=size,Cq1=Cq1,Cq2=Cq2,Nq=Nq,chl=chl,sp=sp))
    end
    for i in N+1:2N
        # agent location
        x = rand(1.5*10:(grid.Nx-0.5)*10)/10
        y = rand(grid.Ny*5:(grid.Ny-0.5)*10)/10
        z = rand(1.5*10:grid.Nz*7.5)/10
        # a normal distribution with mean variance
        radm = max(0.05, rand(Normal(mean,var)))
        gen  = 1
        size = radm
        Cq1  = Cquota[2]*1E-3
        Cq2  = Cquota[2]*radm
        Nq   = 13/106*2*Cq2
        chl  = Cq2*0.4
        sp   = 2
        push!(phyts0,(x=x,y=y,z=z,gen=gen,size=size,Cq1=Cq1,Cq2=Cq2,Nq=Nq,chl=chl,sp=sp))
    end
    B = [phyts0]
    return B
end

# set up initial nutrients fields according to grids
# nut is an array of 6 elements, each element is a kind of nutrient
function setup_nutrients(g,nut)
    DIC = ones(g.Ny, g.Nx, g.Nz)
    DIN = ones(g.Ny, g.Nx, g.Nz)
    DOC = ones(g.Ny, g.Nx, g.Nz)
    DON = ones(g.Ny, g.Nx, g.Nz)
    POC = ones(g.Ny, g.Nx, g.Nz)
    PON = ones(g.Ny, g.Nx, g.Nz)
    DIC .= DIC .* nut[1]
    DIN .= DIN .* nut[2]
    DOC .= DOC .* nut[3]
    DON .= DON .* nut[4]
    POC .= POC .* nut[5]
    PON .= PON .* nut[6]
    nutrients = nutrient_fields(DIC, DIN, DOC, DON, POC, PON)
    return nutrients
end


