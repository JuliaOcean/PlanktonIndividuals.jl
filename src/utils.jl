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
    read_offline_vels(VelOfflineOpt, t)
Read velocity fields of selected grids at 't' from cluster
Return a velocity 'struc'
extra cols & rows are read in for double grids
z starts from the bottom
"""
function read_offline_vels(VelOfflineOpt::Dict,t::Int64)
    vfroot  = VelOfflineOpt["velpath"]
    itList  = VelOfflineOpt["itList"]
    Grid_sel= VelOfflineOpt["GridSel"]
    tN      = VelOfflineOpt["tN"]
    Nx⁻ = Grid_sel["Nx"][1]-1; Nx⁺ = Grid_sel["Nx"][2]+1;
    Ny⁻ = Grid_sel["Ny"][1]-1; Ny⁺ = Grid_sel["Ny"][2]+1;
    Nz⁻ = Grid_sel["Nz"][1]; Nz⁺ = Grid_sel["Nz"][2]+1;
    Nx = Nx⁺ - Nx⁻ + 1; Ny = Ny⁺ - Ny⁻ + 1; Nz = Nz⁺ - Nz⁻ + 2;
    if Nx == 1
        u = zeros(Nx, Ny, Nz)
    else
        fuvel = open(vfroot*"/UVEL/_."*lpad(string(itList[t+tN]),10,"0")*".data")
        uvel = zeros(Float32, 1080, 2700, 90)
        read!(fuvel, uvel); uvel .= ntoh.(uvel)
        close(fuvel);
        u = zeros(Nx, Ny, Nz)
        u[:,:,2:end] = uvel[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺]
        u = reverse(u,dims=3)
    end

    if Ny == 1
        v = zeros(Nx, Ny, Nz)
    else
        fvvel = open(vfroot*"/VVEL/_."*lpad(string(itList[t+tN]),10,"0")*".data")
        vvel = zeros(Float32, 1080, 2700, 90)
        read!(fvvel, vvel); vvel .= ntoh.(vvel)
        close(fvvel);
        v = zeros(Nx, Ny, Nz)
        v[:,:,2:end] = vvel[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺]
        v = reverse(v,dims=3)
    end

    if Nz == 1
        w = zeros(Nx, Ny, Nz+1)
    else
        fwvel = open(vfroot*"/WVEL/_."*lpad(string(itList[t+tN]),10,"0")*".data")
        wvel = zeros(Float32, 1080, 2700, 90)
        read!(fwvel, wvel); wvel .= ntoh.(wvel)
        close(fwvel);
        w = zeros(Nx, Ny, Nz+1)
        w[:,:,3:end] = wvel[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺]
        w = reverse(w,dims=3)
    end
    vel = velocity(u, v, w)
    return vel
end

"""
    grid_offline(GridOfflineOpt)
Read grid information of selected grids from cluster
Return a grid 'struc'
"""
function grid_offline(GridOfflineOpt::Dict)
    fieldroot = GridOfflineOpt["gridpath"]
    Grid_sel  = GridOfflineOpt["GridSel"]
    nx=1080;ny=2700;nz=40;
    fxg = open(fieldroot*"XG.data","r");
    fyg = open(fieldroot*"YG.data","r");
    fxc = open(fieldroot*"XC.data","r");
    fyc = open(fieldroot*"YC.data","r");
    fdx = open(fieldroot*"DXG.data","r");
    fdy = open(fieldroot*"DYG.data","r");
    fdrf= open(fieldroot*"DRF.data","r");
    fdxc= open(fieldroot*"DXC.data","r");
    fdyc= open(fieldroot*"DYC.data","r");
    fdrc= open(fieldroot*"DRC.data","r");
    fAz = open(fieldroot*"RAC.data","r");
    fhfc= open(fieldroot*"hFacC.data","r");
    fhfs= open(fieldroot*"hFacS.data","r");
    fhfw= open(fieldroot*"hFacW.data","r");
    xf = zeros(Float32,nx,ny); yf = zeros(Float32,nx,ny);
    xc = zeros(Float32,nx,ny); yc = zeros(Float32,nx,ny);
    dx = zeros(Float32,nx,ny); dy = zeros(Float32,nx,ny);
    drf= zeros(Float32,nz); drc = zeros(Float32,nz);
    dxc= zeros(Float32,nx,ny); dyc= zeros(Float32,nx,ny);
    Az = zeros(Float32,nx,ny); hFC= zeros(Float32,nx,ny,nz);
    hFS= zeros(Float32,nx,ny,nz);hFW= zeros(Float32,nx,ny,nz);
    read!(fxg,xf); read!(fyg,yf); read!(fxc,xc); read!(fyc,yc);
    read!(fdx,dx); read!(fdy,dy); read!(fdrf,drf); read!(fdxc,dxc);
    read!(fdyc,dyc); read!(fdrc,drc); read!(fAz,Az);
    read!(fhfc,hFC); read!(fhfs,hFS); read!(fhfw,hFW);
    close(fxg);close(fyg);close(fxc);close(fyc);close(fdx);close(fdy);
    close(fdrf);close(fAz);close(fhfc);close(fhfs);close(fhfw);
    xf .= ntoh.(xf); yf .= ntoh.(yf);
    xc .= ntoh.(xc); yc .= ntoh.(yc);
    dx .= ntoh.(dx); dy .= ntoh.(dy);
    drf.= ntoh.(drf);drc.= ntoh.(drc);
    dxc.= ntoh.(dxc);dyc.= ntoh.(dyc);
    Az .= ntoh.(Az); hFC.= ntoh.(hFC);
    hFS.= ntoh.(hFS);hFW.= ntoh.(hFW);
    zf = -cumsum(drf); pushfirst!(zf,0); zc = 0.5*(zf[1:end-1]+zf[2:end]);
    zf = reverse(zf); zc = reverse(zc); drf = reverse(drf); drc = reverse(drc);
    Ax = zeros(nx,ny,nz); Ay = zeros(nx,ny,nz); V = zeros(nx,ny,nz);
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                Ax[i,j,k] = drf[k] * dy[i,j] * hFW[i,j,k]
                Ay[i,j,k] = drf[k] * dx[i,j] * hFS[i,j,k]
                V[i,j,k] = drf[k] * Az[i,j] * hFC[i,j,k]
            end
        end
    end
    # seletc grids
    Nx⁻ = Grid_sel["Nx"][1]; Nx⁺ = Grid_sel["Nx"][2]
    Ny⁻ = Grid_sel["Ny"][1]; Ny⁺ = Grid_sel["Ny"][2]
    Nz⁻ = Grid_sel["Nz"][1]; Nz⁺ = Grid_sel["Nz"][2]
    xcS = xc[Nx⁻:Nx⁺, Ny⁻:Ny⁺]; ycS = yc[Nx⁻:Nx⁺, Ny⁻:Ny⁺];
    xfS = xf[Nx⁻:Nx⁺+1, Ny⁻:Ny⁺]; yfS = yf[Nx⁻:Nx⁺, Ny⁻:Ny⁺+1];
    dxS = dx[Nx⁻:Nx⁺, Ny⁻:Ny⁺]; dyS = dy[Nx⁻:Nx⁺, Ny⁻:Ny⁺];
    dxcS= dxc[Nx⁻:Nx⁺,Ny⁻:Ny⁺];dycS= dyc[Nx⁻:Nx⁺, Ny⁻:Ny⁺];
    zcS = zc[Nz⁻:Nz⁺]; zfS = zf[Nz⁻:Nz⁺+1];
    drfS = drf[Nz⁻:Nz⁺]; drcS = drc[Nz⁻:Nz⁺];
    AzS = Az[Nx⁻:Nx⁺, Ny⁻:Ny⁺]; AxS = Ax[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺];
    AyS = Ay[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺];
    VS  =  V[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺];
    Nx, Ny, Nz = size(VS)
    g = grids(xcS, ycS, zcS, xfS, yfS, zfS, dxS, dyS, drfS, dxcS, dycS, drcS, AxS, AyS, AzS, VS, Nx, Ny, Nz)
    return g
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
    xC = repeat(Ogrid.xC[:,1,1],1,Ny)
    yC = rotl90(repeat(Ogrid.yC[1,:,1],1,Nx))
    zC = collect(Ogrid.zC[1,1,:])
    xF = repeat(Ogrid.xF[:,1,1],1,Ny)
    yF = rotl90(repeat(Ogrid.yF[1,:,1],1,Nx))
    zF = collect(Ogrid.zF[1,1,:])
    dxF =  repeat(collect(Ogrid.Δx),Nx,Ny)
    dyF =  repeat(collect(Ogrid.Δy),Nx,Ny)
    dzF =  repeat(collect(Ogrid.Δz),Nz)
    dxC =  repeat(collect(Ogrid.Δx),Nx,Ny)
    dyC =  repeat(collect(Ogrid.Δy),Nx,Ny)
    dzC =  repeat(collect(Ogrid.Δz),Nz)
    Ax = zeros(Nx, Ny, Nz); Ay = zeros(Nx, Ny, Nz);
    Az = zeros(Nx, Ny); V = zeros(Nx, Ny, Nz)
    for i in 1:Nx
        for j in 1:Ny
            Az[i,j] = dxF[i,j] * dyF[i,j]
            for k in 1:Nz
                Ax[i,j,k] = dzF[k] * dyF[i,j]
                Ay[i,j,k] = dzF[k] * dxF[i,j]
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
    xF = g.xF[:,1]; yF = g.yF[1,:]; zF = g.zF;
    xind = findall(t -> t≤x, xF[1:end-1])[end]
    yind = findall(t -> t≤y, yF[1:end-1])[end]
    zind = findall(t -> t≤z, zF[1:end-1])[end]
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
