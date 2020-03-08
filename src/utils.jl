"""
    read_default_IR_input(nTime, ΔT, grid)
    'nTime' -> number of time step
    'ΔT' -> length of each time step
    'grid' -> grid information
Read input of irradiance from default binary file
"""
function read_IR_input(nTime::Int64, ΔT::Int64,grid,
                       path = dirname(pathof(PhytoAgentModel))*"/../samples/PAR.bin")
    # irradiance(μmol photons/m^2)
    # start from mid-night
    # will change later to make it able to choose different month
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
    PAR_surf = zeros(grid.Nx, grid.Ny, size(PAR,1))
    for i in 1:size(PAR,1)
        PAR_surf[:,:,i] .= PAR[i]
    end
    return PAR_surf
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
                         path = dirname(pathof(PhytoAgentModel))*"/../samples/temp.bin")
    # will change later to make it able to choose different month
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
        temp_domain[:,:,1,i] .= temp[i]
    end
    # vertical temperature gradient
    for j in 2:grid.Nz
        temp_domain[:,:,j,:] .= temp_domain[:,:,j-1,:] .- (∂T∂z*(grid.zC[j-1]-grid.zC[j]))
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
    grid_Ogrids(Ogrid)
Read grid information from Oceananigans
Return a grid 'struc'
'z' starts from the bottom
"""
function read_Ogrids(Ogrid)
    Nx = Ogrid.Nx
    Ny = Ogrid.Ny
    Nz = Ogrid.Nz
    xC = repeat(collect(Ogrid.xC),1,Ny)
    yC = rotl90(repeat(collect(Ogrid.yC),1,Nx))
    zC = collect(Ogrid.zC)
    xF = repeat(collect(Ogrid.xF),1,Ny)
    yF = rotl90(repeat(collect(Ogrid.yF),1,Nx))
    zF = collect(Ogrid.zF)
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
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
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
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
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
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        x,y,z = which_grid(phyt, gird)
        HD[x,y] = HD[x,y] + 1.0
    end
    return HD
end

"""
    write_nut_nc_each_step(g, nut, t)
Write a NetCDF file of nutrient fields at each time step
Default filepath -> "results/nutrients/nut."*lpad(string(t),4,"0")*".nc"
"""
function write_nut_nc_each_step(g::grids, nut::nutrient_fields, t::Int64, filepath::String)
    xC_attr = Dict("longname" => "Locations of the cell centers in the x-direction.", "units" => "m")
    yC_attr = Dict("longname" => "Locations of the cell centers in the y-direction.", "units" => "m")
    zC_attr = Dict("longname" => "Locations of the cell centers in the z-direction.", "units" => "m")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    nccreate(filepath, "DIC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "NH4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "NO3", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "PO4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=P_attr);
    nccreate(filepath, "DOC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "DON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "DOP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=P_attr);
    nccreate(filepath, "POC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "PON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "POP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=P_attr);
    ncwrite(nut.DIC,filepath,"DIC"); ncwrite(nut.NH4,filepath,"NH4");
    ncwrite(nut.NO3,filepath,"NO3"); ncwrite(nut.PO4,filepath,"PO4");
    ncwrite(nut.DOC,filepath,"DOC"); ncwrite(nut.DON,filepath,"DON"); ncwrite(nut.DOP,filepath,"DOP");
    ncwrite(nut.POC,filepath,"POC"); ncwrite(nut.PON,filepath,"PON"); ncwrite(nut.POP,filepath,"POP");
    nothing
end

"""
    write_nut_nc_alltime(a, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime)
Write a NetCDF file of nutrient fields for the whole run, especially for 0D configuration
Default filepath -> "results/nutrients.nc"
"""
function write_nut_nc_alltime(g::grids, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime,
                              filepath = "./results/nutrients.nc")
    tt = collect(1:nTime);
    xC_attr = Dict("longname" => "Locations of the cell centers in the x-direction.", "units" => "m")
    yC_attr = Dict("longname" => "Locations of the cell centers in the y-direction.", "units" => "m")
    zC_attr = Dict("longname" => "Locations of the cell centers in the z-direction.", "units" => "m")
    T_attr = Dict("longname" => "Time", "units" => "H")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    nccreate(filepath, "DIC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "NH4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "NO3", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "PO4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=P_attr);
    nccreate(filepath, "DOC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "DON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "DOP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=P_attr);
    nccreate(filepath, "POC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "PON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "POP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=P_attr);
    ncwrite(DIC,filepath,"DIC"); ncwrite(NH4,filepath,"NH4");
    ncwrite(NO3,filepath,"NO3"); ncwrite(PO4,filepath,"PO4");
    ncwrite(DOC,filepath,"DOC"); ncwrite(DON,filepath,"DON"); ncwrite(DOP,filepath,"DOP");
    ncwrite(POC,filepath,"POC"); ncwrite(PON,filepath,"PON"); ncwrite(POP,filepath,"POP");
    ncclose(filepath)
    return nothing
end

"""
    write_nut_cons(g, gtr, nutₜ, vel, t, filepath)
Compute total gtr (supposed to be 0), and surface vertical tracer flux(supposed to be 0)
Write a brief summary of each time step into a txt file
"""
function write_nut_cons(g::grids, gtr::nutrient_fields, nutₜ::nutrient_fields, vel::velocity, t::Int64, filepath)
    Σgtrⁿ = sum(gtr.NH4 .* g.V)+sum(gtr.NO3 .* g.V)+sum(gtr.DON .* g.V)+sum(gtr.PON .* g.V)
    Σgtrᶜ = sum(gtr.DIC .* g.V)+sum(gtr.DOC .* g.V)+sum(gtr.POC .* g.V)
    Σgtrᵖ = sum(gtr.PO4 .* g.V)+sum(gtr.DOP .* g.V)+sum(gtr.POP .* g.V)
    ΣsurFⁿ= sum((nutₜ.NH4[:,:,1]+nutₜ.NO3[:,:,1]+nutₜ.DON[:,:,1]+nutₜ.PON[:,:,1]) .* g.Az .* vel.w[2:end-1,2:end-1,end-2])
    ΣsurFᶜ= sum((nutₜ.DIC[:,:,1]+nutₜ.DOC[:,:,1]+nutₜ.POC[:,:,1]) .* g.Az .* vel.w[2:end-1,2:end-1,end-2])
    ΣsurFᵖ= sum((nutₜ.PO4[:,:,1]+nutₜ.DOP[:,:,1]+nutₜ.POP[:,:,1]) .* g.Az .* vel.w[2:end-1,2:end-1,end-2])
    Cio = open(filepath*"cons_C.txt","a"); Nio = open(filepath*"cons_N.txt","a");
    Pio = open(filepath*"cons_P.txt","a");
    println(Cio,@sprintf("%3.0f  %.16E  %.16E  %.8E",t,Σgtrᶜ,ΣsurFᶜ,Σgtrᶜ+ΣsurFᶜ))
    println(Nio,@sprintf("%3.0f  %.16E  %.16E  %.8E",t,Σgtrⁿ,ΣsurFⁿ,Σgtrⁿ+ΣsurFⁿ))
    println(Pio,@sprintf("%3.0f  %.16E  %.16E  %.8E",t,Σgtrᵖ,ΣsurFᵖ,Σgtrᵖ+ΣsurFᵖ))
    close(Cio);close(Nio);close(Pio);
end

"""
    write_pop_dynamics(t, agent_num, counts, filepath)
Write a brief summary of population changes at each time step into a txt file
"""
function write_pop_dynamics(t::Int64, pop, counts, filepath)
    POPio = open(filepath*"pop_dynamics.txt","a");
    println(POPio,@sprintf("%3.0f  %7.0f  %5.0f  %5.0f  %5.0f",
                           t,pop,counts.divid, counts.graze,counts.death))
    close(POPio);
end

"""
    write_output(individuals,filepath,time)
write model output of individuals at each time step in a binary file
time = model.t*ΔT
"""
function write_output(individuals, filepath, time)
    phytos = individuals.phytos
    path = filepath*"phy_"*lpad(time, 10, "0")*".bin"
    if individuals.zoos == nothing
        open(path, "w") do io
            serialize(io, phytos)
        end
    else
        open(path, "w") do io
            serialize(io, phytos)
        end
        path_zoo = filepath*"zoo_"*lpad(time, 10, "0")*".bin"
        open(path_zoo, "w") do io
            serialize(io, individuals.zoos)
        end
    end
end

"""
    PrepRunDir(res::String="results/")
Create `res/` folder if needed. Remove old files from it if needed.
"""
function PrepRunDir(res::String="./results/")
    isdir(res) && rm(res, recursive=true)
    mkdir(res)
    mkdir("$res"*"nutrients/")
    return res
end
