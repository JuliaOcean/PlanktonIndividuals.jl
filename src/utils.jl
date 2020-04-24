"""
    read_default_IR_input(nTime, ΔT, grid)
    'nTime' -> number of time step
    'ΔT' -> length of each time step
    'grid' -> grid information
Read input of irradiance from default binary file
"""
function read_IR_input(nTime::Int64, ΔT::Int64,grid,
                       path = dirname(pathof(PlanktonIndividuals))*"/../samples/PAR.bin")
    # irradiance(μmol photons/s/m^2)
    # start from mid-night
    PAR_hour = deserialize(path)
    # convert hour to second in a day
    t_htos = collect(0:1:24) .* 3600
    # convert each time step to second in a day
    t_ΔT = collect(1:ΔT:86400)
    # interpolation of PAR time series
    itp_PAR = interpolate((t_htos,), PAR_hour, Gridded(Linear()));
    PAR_itped = itp_PAR.(t_ΔT)
    # deal with nTime
    if nTime*ΔT < 86400
        PAR = PAR_itped[1:nTime]
    else
        nday = (nTime*ΔT)÷86400
        res = (nTime*ΔT)%86400÷ΔT
        PAR = repeat(PAR_itped,nday)
        PAR = vcat(PAR,PAR_itped[1:res])
    end
    # expand to the whole domain surface
    PAR_domain = zeros(grid.Nx, grid.Ny, grid.Nz, size(PAR,1))
    for i in 1:size(PAR,1)
        PAR_domain[:,:,end,i] .= PAR[i]
    end
    return PAR_domain
end

"""
    read_default_temp_input(nTime, ΔT, grid, ∂T∂z)
    'nTime' -> number of time step
    'ΔT' -> length of each time step
    'grid' -> grid information
    '∂T∂z' -> linear vertical temp gradient
Read input of temperature from default binary file
"""
function read_temp_input(nTime::Int64, ΔT::Int64, grid, ∂T∂z=0.04,
                         path = dirname(pathof(PlanktonIndividuals))*"/../samples/temp.bin")
    temp_hour = deserialize(path)
    # convert hour to second in a day
    t_htos = collect(0:1:24) .* 3600
    # convert each time step to second in a day
    t_ΔT = collect(1:ΔT:86400)
    # interpolation of PAR time series
    itp_temp = interpolate((t_htos,), temp_hour, Gridded(Linear()));
    temp_itped = itp_temp.(t_ΔT)
    # deal with nTime
    if nTime*ΔT < 86400
        temp = temp_itped[1:nTime]
    else
        nday = (nTime*ΔT)÷86400
        res = (nTime*ΔT)%86400÷ΔT
        temp = repeat(temp_itped,nday)
        temp = vcat(temp,temp_itped[1:res])
    end
    # espand to the whole domain
    temp_domain = zeros(grid.Nx, grid.Ny, grid.Nz, size(temp,1))
    for i in 1:size(temp,1)
        temp_domain[:,:,end,i] .= temp[i]
    end
    # vertical temperature gradient
    for j in grid.Nz-1:-1:1
        temp_domain[:,:,j,:] .= temp_domain[:,:,j+1,:] .- (∂T∂z*(grid.zC[j+1]-grid.zC[j]))
    end
    return temp_domain
end


"""
    update_params!(parameters, tmp)
Update parameter values based on .yaml file
'parameters' is default parameter set
'tmp' is the parameters need to update
"""
function update_params!(parameters::Dict, tmp::Dict)
    tmp_keys = collect(keys(tmp))
    for key in tmp_keys
        if length(findall(x->x==key, collect(keys(parameters))))==0
            print("PARAM: parameter not found \n")
        else
            parameters[key] = tmp[key]
        end
    end
    return parameters
end

"""
    grid_Ogrids(Ogrid)
Read grid information from Oceananigans
Return a grid 'struc'
'z' starts from the bottom
"""
function read_Ogrids(Ogrid)
    Nx = Ogrid.Nx
    Ny = Ogrid.Ny
    Nz = Ogrid.Nz
    xC = collect(Ogrid.xC[:,1,1])
    yC = collect(Ogrid.yC[1,:,1])
    zC = collect(Ogrid.zC[1,1,:])
    xF = collect(Ogrid.xF[:,1,1])
    yF = collect(Ogrid.yF[1,:,1])
    zF = collect(Ogrid.zF[1,1,:])
    dxF =  repeat(collect(Ogrid.Δx),Nx)
    dyF =  repeat(collect(Ogrid.Δy),Ny)
    dzF =  repeat(collect(Ogrid.Δz),Nz)
    dxC =  repeat(collect(Ogrid.Δx),Nx)
    dyC =  repeat(collect(Ogrid.Δy),Ny)
    dzC =  repeat(collect(Ogrid.Δz),Nz)
    Ax = zeros(Nx, Ny, Nz); Ay = zeros(Nx, Ny, Nz);
    Az = zeros(Nx, Ny); V = zeros(Nx, Ny, Nz)
    for i in 1:Nx
        for j in 1:Ny
            Az[i,j] = dxF[i] * dyF[j]
            for k in 1:Nz
                Ax[i,j,k] = dzF[k] * dyF[j]
                Ay[i,j,k] = dzF[k] * dxF[i]
                V[i,j,k] = dzF[k] * Az[i,j]
            end
        end
    end
    g = grids(xC, yC, zC, xF, yF, zF, dxF, dyF, dzF, dxC, dyC, dzC, Ax, Ay, Az, V, Nx, Ny, Nz)
    return g
end

"""
    which_grid(phyt,g)
decide which grid an individual is in
return grid indices
"""
function which_grid(phyt, g)
    x = phyt[1]; y = phyt[2]; z = phyt[3];
    xind = findall(t -> t≤x, g.xF[2:end-1])[end]
    yind = findall(t -> t≤y, g.yF[2:end-1])[end]
    zind = findall(t -> t≤z, g.zF[2:end-1])[end]
    return xind, yind, zind
end

"""
    count_chl(phyts_a, grid)
Compute total chl concentration in each grid cell
"""
function count_chl(phyts_a, grid)
    cells = zeros(grid.Nx, grid.Ny, grid.Nz)
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        x,y,z = which_grid(phyt, grid)
        cells[x, y, z] = cells[x, y, z] + phyt[12]
    end
    cells .= cells ./ grid.V
    return cells
end

"""
    count_vertical_num(phyts_a)
Count individual numbers in each meter vertically, accumulate horizontally
"""
function count_vertical_num(phyts_a)
    VD = zeros(500)
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        z = -trunc(Int, phyt[3]) + 1
        VD[z] = VD[z] + 1.0
    end
    return VD
end

"""
    count_horizontal_num(phyts_a, grid)
Count individual numbers in each horizontal grid cell, accumulate in vertical direction
"""
function count_horizontal_num(phyts_a,grid)
    HD = zeros(grid.Nx,grid.Ny)
    for i in 1:size(phyts_a,2)
        phyt = phyts_a[:,i]
        x,y,z = which_grid(phyt, gird)
        HD[x,y] = HD[x,y] + 1.0
    end
    return HD
end
