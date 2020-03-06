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
function find_feeding_area(zplk::Array, phyts::Array, radius::Float64)
    if size(phyts,1) == 0
        return (false, nothing, (0.0,0.0,0.0))
    end
    dist= (phyts[1] .- zplk[1]).^2 .+ (phyts[2] .- zplk[2]).^2 .+ (phyts[3] .- zplk[3]).^2
    # dist_df = DataFrame(ID=collect(1:1:size(phyts,1)), dist=dist)
    dist_df = [collect(1:1:size(phyts,1)) dist]
    dist_df = sortslices(dist_df, dims = 1, by = x -> x[2])
    indices = findall(x -> x<radius, dist_df[:,2])
    if size(indices,1) == 0
        if dist_df[1,2] > (radius * 1.0e6)
            return (false, dist_df[1,1], nothing)
        else
            x = copy(phyts[dist_df[1,1],1])
            y = copy(phyts[dist_df[1,1],2])
            z = copy(phyts[dist_df[1,1],3])
            return (false, dist_df[1,1], (x,y,z))
        end
    elseif size(indices,1) == size(dist_df,1)
        return (true, dist_df[indices, 1], nothing)
    else
        nearest_ID = dist_df[1,indices[end]+1]
        x = copy(phyts[nearest_ID,1])
        y = copy(phyts[nearest_ID,2])
        z = copy(phyts[nearest_ID,3])
        return (true, dist_df[indices, 1], (x,y,z))
    end
end

"""
    chase_prey(zplk, cord, travel_dist, grid)
Zooplankton chase after the nearest prey
return the distance the predator traveled
"""
function chase_prey(zplk::Array, cord, travel_dist::Float64, grid)
    dx = cord[1] - zplk[1]
    dy = cord[2] - zplk[2]
    dz = cord[3] - zplk[3]
    dist = dx^2 + dy^2 + dz^2
    dratio = travel_dist/dist
    if dratio ≥ 1
        zplk[1] = zplk[1] + dx
        zplk[2] = zplk[2] + dy
        zplk[3] = zplk[3] + dz
        return dist
    else
        zplk[1] = zplk[1] + dx * sqrt(dratio)
        zplk[2] = zplk[2] + dy * sqrt(dratio)
        zplk[3] = zplk[3] + dz * sqrt(dratio)
        # periodic domain
        zplk[3] = max(grid.zF[1],min(grid.zF[end],zplk[3] ))
        zplk[1] = periodic_domain(grid.xF, zplk[1])
        zplk[2] = periodic_domain(grid.yF, zplk[2])
        return travel_dist
    end
end

"""
    rand_walk(zplk::DataFrameRow, grid, travel_dist)
If zooplankton cannot find any prey, it will swim randomly
"""
function rand_walk!(zplk::Array, grid, travel_dist)
    dx = rand(Uniform(-1,1))
    dy = rand(Uniform(-1,1))
    dz = rand(Uniform(-1,1))
    zplk[1] += dx*travel_dist
    zplk[2] += dy*travel_dist
    zplk[3] += dz*travel_dist
    # periodic domain
    zplk[3] = max(grid.zF[1],min(grid.zF[end],zplk[3] ))
    zplk[1] = periodic_domain(grid.xF, zplk[1])
    zplk[2] = periodic_domain(grid.yF, zplk[2])
    return nothing
end


"""
    sloppy_feed(zplk, phyts_feed, params)
zooplankton feed on selected phytoplankton
return C and N content as exports to environment
"""
function sloppy_feed(zplk::Array, phyts_feed::Array, params)
    TOC = sum(phyts_feed[:,8]) + sum(phyts_feed[:,9])
    TON = sum(phyts_feed[:,10])
    TOP = sum(phyts_feed[:,11])
    Feed_CN = (params["slpyFracC"] * TOC) / (params["slpyFracN"] * TON)
    Feed_CP = (params["slpyFracC"] * TOC) / (params["slpyFracP"] * TOP)
    Zoo_CN = zplk[9] / zplk[10]
    Zoo_CP = zplk[9] / zplk[11]
    if Feed_CN ≥ Zoo_CN
        dCq2 = params["slpyFracN"] * TON * Zoo_CN
        dsize= dCq2 / zplk[9]
    else
        dCq2 = params["slpyFracC"] * TOC
        dsize= dCq2 / zplk[9]
    end
    if Feed_CP ≥ Zoo_CP
        dCq2 = params["slpyFracP"] * TOP * Zoo_CP
        dsize= dCq2 / zplk[9]
    else
        dCq2 = params["slpyFracC"] * TOC
        dsize= dCq2 / zplk[9]
    end
    zplk[9]  = zplk[9] + dCq2
    zplk[10] = zplk[10]  + dCq2 / Zoo_CN
    zplk[11] = zplk[11]  + dCq2 / Zoo_CP
    zplk[7]  = zplk[7]+ dsize
    Cexport = TOC - dCq2
    Nexport = TON - dCq2 / Zoo_CN
    Pexport = TOP - dCq2 / Zoo_CP
    return Cexport, Nexport, Pexport
end


"""
    zoo_update(zoos, phyts, ΔT)
update individuals of zooplankton one time step forward
both physiology processes and active movement
Control flow: Death -> Graze & Swim & Grow -> Reproduction
"""
function zoo_update(model, ΔT::Int64)
    g = model.grid
    params = model.params
    temp = model.temp
    phyts = model.individuals.phytos
    zoos = model.individuals.zoos
    counts = pop_counts()
    #set up a dataframe to record all updated agents
    zoos_b = zeros(1,size(zoos,2))

    consume = nutrients_init(g)
    for i in 1:size(zoos,1)
        zplk = zoos[i,:]
        x, y, z = which_grid(zplk, g)

        # compute death probability after a certain age
        reg_age = max(0.0, zplk[6] - params["death_age"])
        shape_factor_death = params["a_death"]*reg_age^params["b_death"]
        P_death = rand(Bernoulli(shape_factor_death/(1+shape_factor_death)))

        if P_death == false
            radius = zplk[7] * 1.5 * 5.0e-5
            feeding = find_feeding_area(zplk, phyts, radius)
            # the longest distance predator can swim
            travel_dist = params["v_zoo"] * ΔT

            # find prey or not
            if feeding[1] == false
                if feeding[3] == nothing
                    rand_walk!(zplk, g, travel_dist)
                else
                    travel_dist = chase_prey(zplk, feeding[3], travel_dist, g)
                end
            elseif feeding[1] == true
                phyts_feed = phyts[feeding[2],:]
                phyts = phyts[setdiff(1:end,feeding[2]),:]
                Cexport, Nexport, Pexport = sloppy_feed(zplk, phyts_feed, params)
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + Cexport * params["grazFracC"]
                consume.DON[x, y, z] = consume.DON[x, y, z] + Nexport * params["grazFracN"]
                consume.DOP[x, y, z] = consume.DOP[x, y, z] + Pexport * params["grazFracP"]
                consume.POC[x, y, z] = consume.POC[x, y, z] + Cexport * (1.0 - params["grazFracC"])
                consume.PON[x, y, z] = consume.PON[x, y, z] + Nexport * (1.0 - params["grazFracN"])
                consume.POP[x, y, z] = consume.POP[x, y, z] + Pexport * (1.0 - params["grazFracP"])
                counts.graze = counts.graze + size(phyts_feed,1)
                graz_ct = graz_ct + size(phyts_feed,1)
                if feeding[3] == nothing
                    rand_walk!(zplk, g, travel_dist)
                else
                    travel_dist = chase_prey(zplk, feeding[3], travel_dist, g)
                end
            end
            zplk[6] = zplk[6] + 1.0*(ΔT/3600)

            # metabolic cost, respiration & swim cost
            mtcost = travel_dist * 0.01 * zplk[9]
            zplk[9] = zplk[9] - mtcost

            # compute probabilities of reproduction

            push!(zoos_b,zplk)
        else
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + zplk[9] * params["mortFracC"]
            consume.DON[x, y, z] = consume.DON[x, y, z] + zplk[10]  * params["mortFracN"]
            consume.DOP[x, y, z] = consume.DOP[x, y, z] + zplk[11]  * params["mortFracP"]
            consume.POC[x, y, z] = consume.POC[x, y, z] + zplk[9] * (1.0 - params["mortFracC"])
            consume.PON[x, y, z] = consume.PON[x, y, z] + zplk[10]  * (1.0 - params["mortFracN"])
            consume.POP[x, y, z] = consume.POP[x, y, z] + zplk[11]  * (1.0 - params["mortFracP"])
        end
    end
    return zoos_b[2:end,:], counts, consume
end
