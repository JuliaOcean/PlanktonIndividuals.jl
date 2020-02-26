###################################################################################
# Grazing by zooplankton, sloppy feeding, and physiology processes of zooplankton #
###################################################################################
"""
    find_feeding_area(zplk, phyts, radius)
'zplk' is a individual of zooplankton
'phyts' is the collection of all phytoplankton at current time step
'radius' is the radius of a sphere represent the feeding area
The function finds all preys within the sphere
Return 'false' if no prey is in this sphere with the nearest prey
Return 'true' if there is at least on prey in this sphere along with
a dataframe containing all preys (temporal ID and rd) and the nearest prey
"""
function find_feeding_area(zplk::DataFrameRow, phyts::DataFrame, radius::Float64)
    if size(phyts,1) == 0
        return (false, nothing, (0.0,0.0,0.0))
    end
    dist= (phyts.x .- zplk.x).^2 .+ (phyts.y .- zplk.y).^2 .+ (phyts.z .- zplk.z).^2
    dist_df = DataFrame(ID=collect(1:1:size(phyts,1)), dist=dist)
    dist_df = sort(dist_df,:dist)
    indices = findall(x -> x<radius, dist_df.dist)
    if size(indices,1) == 0
        x = copy(phyts[dist_df.ID[1],:].x)
        y = copy(phyts[dist_df.ID[1],:].y)
        z = copy(phyts[dist_df.ID[1],:].z)
        return (false, dist_df[1,:], (x,y,z))
    elseif size(indices,1) == size(dist_df,1)
        return (true, dist_df[indices, :], (0.0,0.0,0.0))
    else
        nearest_ID = dist_df.ID[indices[end]+1]
        x = copy(phyts[nearest_ID,:].x)
        y = copy(phyts[nearest_ID,:].y)
        z = copy(phyts[nearest_ID,:].z)
        return (true, dist_df[indices, :], (x,y,z))
    end
end

"""
    chase_prey(zplk, cord, travel_dist, grid)
Zooplankton chase after the nearest prey
return the distance the predator traveled
"""
function chase_prey(zplk::DataFrameRow, cord, travel_dist::Float64, grid)
    dx = cord[1] - zplk.x
    dy = cord[2] - zplk.y
    dz = cord[3] - zplk.z
    dist = dx^2 + dy^2 + dz^2
    dratio = travel_dist/dist
    if dratio ≥ 1
        zplk.x = zplk.x + dx
        zplk.y = zplk.y + dy
        zplk.z = zplk.z + dz
        return dist
    else
        zplk.x = zplk.x + dx * sqrt(dratio)
        zplk.y = zplk.y + dy * sqrt(dratio)
        zplk.z = zplk.z + dz * sqrt(dratio)
        # periodic domain
        zplk.z = max(grid.zF[end],min(grid.zF[1],zplk.z ))
        zplk.x = periodic_domain(grid.xF, zplk.x)
        zplk.y = periodic_domain(grid.yF, zplk.y)
        return travel_dist
    end
end


"""
    sloppy_feed(zplk, phyts_feed, params)
zooplankton feed on selected phytoplankton
return C and N content as exports to environment
"""
function sloppy_feed(zplk::DataFrameRow, phyts_feed::DataFrame, params)
    TOC = sum(phyts_feed.Cq1) + sum(phyts_feed.Cq2)
    TON = sum(phyts_feed.Nq)
    Feed_CN = (params["slpyFracC"] * TOC) / (params["slpyFracN"] * TON)
    Zoo_CN = zplk.Cq2 / zplk.Nq
    if Feed_CN ≥ Zoo_CN
        dCq2 = params["slpyFracN"] * TON * Zoo_CN
        dsize= dCq2 / zplk.Cq2
    else
        dCq2 = params["slpyFracC"] * TOC
        dsize= dCq2 / zplk.Cq2
    end
    zplk.Cq2 = zplk.Cq2 + dCq2
    zplk.Nq  = zplk.Nq  + dCq2 / Zoo_CN
    zplk.size= zplk.size+ dsize
    Cexport = TOC - dCq2
    Nexport = TON - dCq2 / Zoo_CN
    return Cexport, Nexport
end


"""
    zoo_update(zoos, phyts, ΔT)
update individuals of zooplankton one time step forward
both physiology processes and active movement
Control flow: Death -> Graze & Swim & Grow -> Reproduction
"""
function zoo_update(zoos::DataFrame, phyts::DataFrame, ΔT::Int64, model)
    g = model.grid
    params = model.params
    temp = model.temp
    dvid_ct = 0; graz_ct = 0; death_ct = 0
    #set up a dataframe to record all updated agents
    zoos_b = DataFrame(x=Float64[], y=Float64[], z=Float64[],
                      gen=Int64[], size=Float64[], Cq1=Float64[],
                      Cq2=Float64[], Nq=Float64[], chl=Float64[],
                      sp=Int64[], age=Float64[])

    consume = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                              zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                              zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    for i in 1:size(zoos,1)
        zplk = zoos[i,:]
        x, y, z = which_grid(zplk, g)

        # compute death probability after a certain age
        reg_age = max(0.0, zplk.age - params["death_age"])
        shape_factor_death = params["a_death"]*reg_age^params["b_death"]
        P_death = rand(Bernoulli(shape_factor_death/(1+shape_factor_death)))

        if P_death == false
            radius = zplk.size * 1.5 * 5.0e-5
            feeding = find_feeding_area(zplk, phyts, radius)
            # the longest distance predator can swim
            travel_dist = params["v_zoo"] * ΔT

            # find prey or not
            if feeding[1] == false
                travel_dist = chase_prey(zplk, feeding[3], travel_dist, g)
            elseif feeding[1] == true
                phyts_feed = phyts[feeding[2].ID,:]
                deleterows!(phyts,sort(feeding[2].ID))
                Cexport,Nexport = sloppy_feed(zplk, phyts_feed, params)
                travel_dist = chase_prey(zplk, feeding[3], travel_dist, g)
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + Cexport * params["grazFracC"]
                consume.DON[x, y, z] = consume.DON[x, y, z] + Nexport * params["grazFracN"]
                consume.POC[x, y, z] = consume.POC[x, y, z] + Cexport * (1.0 - params["grazFracC"])
                consume.PON[x, y, z] = consume.PON[x, y, z] + Nexport * (1.0 - params["grazFracN"])
                graz_ct = graz_ct + size(phyts_feed,1)
            end
            zplk.age = zplk.age + 1.0*(ΔT/3600)

            # metabolic cost, respiration & swim cost
            mtcost = travel_dist * 0.01 * zplk.Cq2
            zplk.Cq2 = zplk.Cq2 - mtcost

            # compute probabilities of reproduction

            push!(zoos_b,zplk)
        else
            death_ct += 1
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + zplk.Cq2 * params["mortFracC"]
            consume.DON[x, y, z] = consume.DON[x, y, z] + zplk.Nq  * params["mortFracN"]
            consume.POC[x, y, z] = consume.POC[x, y, z] + zplk.Cq2 * (1.0 - params["mortFracC"])
            consume.PON[x, y, z] = consume.PON[x, y, z] + zplk.Nq  * (1.0 - params["mortFracN"])
        end
    end
    return zoos_b, (dvid_ct, graz_ct, death_ct), consume
end
