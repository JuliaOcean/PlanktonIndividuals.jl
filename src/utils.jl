"""
    read_default_IR_input(nTime, ΔT, grid)
    nTime -> number of time step
    ΔT -> length of each time step
    grid -> grid information
Read input of irradiance from default binary file
"""
function read_default_IR_input(nTime::Int64, ΔT::Int64,grid)
    # irradiance(μmol photons/m^2)
    # start from mid-night
    # will change later to make it able to choose different month
    path = dirname(pathof(PhytoAgentModel))*"/../samples/PAR.bin"
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
    nTime -> number of time step
    ΔT -> length of each time step
    grid -> grid information
    ∂T∂z -> linear vertical temp gradient
Read input of temperature from default binary file
"""
function read_default_temp_input(nTime::Int64, ΔT::Int64, grid, ∂T∂z=0.04)
    # will change later to make it able to choose different month
    path = dirname(pathof(PhytoAgentModel))*"/../samples/temp.bin"
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
"""
function read_offline_vels(VelOfflineOpt::Dict,t::Int64)
    vfroot  = VelOfflineOpt["velpath"]
    itList  = VelOfflineOpt["itList"]
    Grid_sel= VelOfflineOpt["GridSel"]
    tN      = VelOfflineOpt["tN"]
    Nx⁻ = Grid_sel["Nx"][1]; Nx⁺ = Grid_sel["Nx"][2]
    Ny⁻ = Grid_sel["Ny"][1]; Ny⁺ = Grid_sel["Ny"][2]
    Nz⁻ = Grid_sel["Nz"][1]; Nz⁺ = Grid_sel["Nz"][2]
    Nx = Nx⁺ - Nx⁻ + 1; Ny = Ny⁺ - Ny⁻ + 1; Nz = Nz⁺ - Nz⁻ + 1;
    if Nx == 1 
        u = zeros(Nx, Ny, Nz)
    else
        fuvel = open(vfroot*"/UVEL/_."*lpad(string(itList[t+tN]),10,"0")*".data")
        uvel = zeros(Float32, 1080, 2700, 90)
        read!(fuvel, uvel); uvel .= ntoh.(uvel)
        close(fuvel);
        u = uvel[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺]
    end

    if Ny == 1
        v = zeros(Nx, Ny, Nz)
    else
        fvvel = open(vfroot*"/VVEL/_."*lpad(string(itList[t+tN]),10,"0")*".data")
        vvel = zeros(Float32, 1080, 2700, 90)
        read!(fvvel, vvel); vvel .= ntoh.(vvel)
        close(fvvel);
        v = vvel[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺]
    end

    if Nz == 1
        w = zeros(Nx, Ny, Nz)
    else
        fwvel = open(vfroot*"/WVEL/_."*lpad(string(itList[t+tN]),10,"0")*".data")
        wvel = zeros(Float32, 1080, 2700, 90)
        read!(fwvel, wvel); wvel .= ntoh.(wvel)
        close(fwvel);
        w = wvel[Nx⁻:Nx⁺, Ny⁻:Ny⁺, Nz⁻:Nz⁺]
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
    xfS = xf[Nx⁻:Nx⁺, Ny⁻:Ny⁺]; yfS = yf[Nx⁻:Nx⁺, Ny⁻:Ny⁺];
    dxS = dx[Nx⁻:Nx⁺, Ny⁻:Ny⁺]; dyS = dy[Nx⁻:Nx⁺, Ny⁻:Ny⁺];
    dxcS= dxc[Nx⁻:Nx⁺,Ny⁻:Ny⁺];dycS= dyc[Nx⁻:Nx⁺, Ny⁻:Ny⁺];
    zcS = zc[Nz⁻:Nz⁺]; zfS = zf[Nz⁻:Nz⁺];
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
"""
function read_Ogrids(Ogrid)
    Nx = Ogrid.Nx
    Ny = Ogrid.Ny
    Nz = Ogrid.Nz
    xC = repeat(collect(Ogrid.xC),1,Ny)
    yC = rotl90(repeat(collect(Ogrid.yC),1,Nx))
    zC = reverse(collect(Ogrid.zC))
    xF = repeat(collect(Ogrid.xF)[1:end-1],1,Ny)
    yF = rotl90(repeat(collect(Ogrid.yF)[1:end-1],1,Nx))
    zF = reverse(collect(Ogrid.zF)[1:end-1])
    Lx =  repeat(collect(Ogrid.Δx),Nx,Ny)
    Ly =  repeat(collect(Ogrid.Δy),Nx,Ny)
    Lz =  repeat(collect(Ogrid.Δz),Nz)
    dxC =  repeat(collect(Ogrid.Δx),Nx,Ny)
    dyC =  repeat(collect(Ogrid.Δy),Nx,Ny)
    dzC =  repeat(collect(Ogrid.Δz),Nz)
    Ax = zeros(Nx, Ny, Nz); Ay = zeros(Nx, Ny, Nz);
    Az = zeros(Nx, Ny); V = zeros(Nx, Ny, Nz)
    for i in 1:Nx
        for j in 1:Ny
            Az[i,j] = Lx[i,j] * Ly[i,j]
            for k in 1:Nz
                Ax[i,j,k] = Lz[k] * Ly[i,j]
                Ay[i,j,k] = Lz[k] * Lx[i,j]
                V[i,j,k] = Lz[k] * Az[i,j]
            end
        end
    end
    g = grids(xC, yC, zC, xF, yF, zF, Lx, Ly, Lz, dxC, dyC, dzC, Ax, Ay, Az, V, Nx, Ny, Nz)
    return g
end

"""
    creat_output(B)
Creat a dataframe to record the average attributes of indivuduals of each time step
"""
function create_output(B::Array{DataFrame,1})
    output = DataFrame(time=0, gen_ave=mean(B[1].gen), spec_ave = mean(B[1].sp),
                       Cq1_ave=mean(B[1].Cq1), Cq2_ave=mean(B[1].Cq2),
                       Nq_ave=mean(B[1].Nq), size_ave=mean(B[1].size),
                       chl_ave=mean(B[1].chl), Population=size(B[1],1),
                       age_ave=mean(B[1].age), dvid=0, graz=0,death=0)
    return output
end

"""
    write_output(t, phyts_b, dvid_ct, graz_ct, death_ct, output)
Compute the average attributes of individuals at 't' time step and push them into 'output' dataframe
"""
function write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,output)
    # summary of current step
    gen_ave=mean(phyts_b.gen)
    spec_ave=mean(phyts_b.sp)
    Cq1_ave=mean(phyts_b.Cq1)
    Cq2_ave=mean(phyts_b.Cq2)
    Nq_ave=mean(phyts_b.Nq)
    size_ave=mean(phyts_b.size)
    chl_ave=mean(phyts_b.chl)
    age_ave=mean(phyts_b.age)
    push!(output,(time=t, gen_ave=gen_ave, spec_ave=spec_ave,
                  Cq1_ave=Cq1_ave, Cq2_ave=Cq2_ave, Nq_ave=Nq_ave,
                  size_ave=size_ave, chl_ave=chl_ave, Population=size(phyts_b,1),
                  age_ave=age_ave, dvid=dvid_ct, graz=graz_ct, death=death_ct))
    return output
end

"""
    count_chl(phyts_a, grid)
Compute total chl concentration in each grid cell
"""
function count_chl(phyts_a, grid)
    cells = zeros(grid.Nx, grid.Ny, grid.Nz)
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        x = trunc(Int, phyt.x)
        y = trunc(Int, phyt.y)
        z = trunc(Int, phyt.z)
        cells[x, y, z] = cells[x, y, z] + phyt.chl
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
        z = -trunc(Int, phyt.z) + 1
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
        x = trunc(Int, phyt.x)
        y = trunc(Int, phyt.y)
        HD[x,y] = HD[x,y] + 1.0
    end
    return HD
end

"""
    convert_coordinates(phyts, grid)
Convert grid indices of each individual into lat, lon, and depth
"""
function convert_coordinates(phyts, grid)
    for i in 1:size(phyts,1)
    phyt = phyts[i,:]
    z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
    dz = phyt.z - z; dx = phyt.x - x; dy = phyt.y - y;
    phyt.x = grid.xF[x,1] + dx * grid.Δx[x,y];
    phyt.y = grid.yF[1,y] + dy * grid.Δy[x,y];
    phyt.z = grid.zF[z] - dz * grid.Lz[z];
    end
end

"""
    sort_species(Bi, Bsp, sp)
Sort individuals of different species into different dataframes
'Bi' is a dataframe row, 'Nsp' is an empty array
"""
function sort_species(Bi, Bsp, sp::Int64)
    phyts = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[],
                      size=Float64[], Cq1=Float64[], Cq2=Float64[],
                      Nq=Float64[], chl=Float64[],sp=Int64[],age=Float64[])
    for j in 1:size(Bi,1)
        if Bi[j,:].sp == sp
            push!(phyts,Bi[j,:])
        end
    end
    push!(Bsp,phyts)
    return Bsp
end

"""
    compute_mean_species(Bsp, nTime)
Compute the average attributes of individuals of different species
"""
function compute_mean_species(Bsp, nTime)
    output = DataFrame(time=Int64[], gen_ave=Float64[], Cq1_ave=Float64[],
                       Cq2_ave=Float64[], Nq_ave=Float64[], size_ave=Float64[],
                       chl_ave=Float64[], Population=Int64[], age_ave=Float64[]);
    for i in 1:nTime
        gen_ave1=mean(Bsp[i].gen)
        Cq1_ave1=mean(Bsp[i].Cq1)
        Cq2_ave1=mean(Bsp[i].Cq2)
        Nq_ave1=mean(Bsp[i].Nq)
        size_ave1=mean(Bsp[i].size)
        chl_ave1=mean(Bsp[i].chl)
        age_ave1=mean(Bsp[i].age)
        push!(output,(time=i, gen_ave=gen_ave1, Cq1_ave=Cq1_ave1,
                      Cq2_ave=Cq2_ave1, Nq_ave=Nq_ave1, size_ave=size_ave1,
                      chl_ave=chl_ave1, Population=size(Bsp[i],1),age_ave=age_ave1))
    end
    return output
end

"""
    write_nut_nc_each_step(g, nut, t)
Write a NetCDF file of nutrient fields at each time step
"""
function write_nut_nc_each_step(g::grids, nut::nutrient_fields, t::Int64)
    filepath = "results/nutrients/nut."*lpad(string(t),4,"0")*".nc"
    xC_attr = Dict("longname" => "Locations of the cell centers in the x-direction.", "units" => "m")
    yC_attr = Dict("longname" => "Locations of the cell centers in the y-direction.", "units" => "m")
    zC_attr = Dict("longname" => "Locations of the cell centers in the z-direction.", "units" => "m")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    isfile(filepath) && rm(filepath)
    nccreate(filepath, "DIC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "DIN", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "DOC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "DON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "POC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "PON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    ncwrite(nut.DIC,filepath,"DIC"); ncwrite(nut.DIN,filepath,"DIN");
    ncwrite(nut.DOC,filepath,"DOC"); ncwrite(nut.DON,filepath,"DON");
    ncwrite(nut.POC,filepath,"POC"); ncwrite(nut.PON,filepath,"PON");
    ncclose(filepath)
    return nothing
end

"""
    write_nut_nc_alltime(a, DIC, DIN, DOC, DON, POC, PON, nTime)
Write a NetCDF file of nutrient fields for the whole run, especially for 0D configuration
"""
function write_nut_nc_alltime(g::grids, DIC, DIN, DOC, DON, POC, PON, nTime)
    filepath = "results/nutrients.nc"
    tt = collect(1:nTime);
    xC_attr = Dict("longname" => "Locations of the cell centers in the x-direction.", "units" => "m")
    yC_attr = Dict("longname" => "Locations of the cell centers in the y-direction.", "units" => "m")
    zC_attr = Dict("longname" => "Locations of the cell centers in the z-direction.", "units" => "m")
    T_attr = Dict("longname" => "Time", "units" => "H")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    isfile(filepath) && rm(filepath)
    nccreate(filepath, "DIC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "DIN", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "DOC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "DON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "POC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "PON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    ncwrite(nut.DIC,filepath,"DIC"); ncwrite(nut.DIN,filepath,"DIN");
    ncwrite(nut.DOC,filepath,"DOC"); ncwrite(nut.DON,filepath,"DON");
    ncwrite(nut.POC,filepath,"POC"); ncwrite(nut.PON,filepath,"PON");
    ncclose(filepath)
    return nothing
end

"""
    write_nut_cons(g, gtr, nutₜ, vel, agent_num, t, death_ct, graz_ct, dvid_ct)
Compute total gtr (supposed to be 0), and surface vertical tracer flux(supposed to be 0)
Write a brief summary of each time step into a txt file
"""
function write_nut_cons(g::grids, gtr::nutrient_fields, nutₜ::nutrient_fields, vel::velocity, agent_num::Int64, t::Int64, death_ct::Int64, graz_ct::Int64, dvid_ct::Int64)
    Σgtrⁿ = sum(gtr.DIN .* g.V)+sum(gtr.DON .* g.V)+sum(gtr.PON .* g.V)
    Σgtrᶜ = sum(gtr.DIC .* g.V)+sum(gtr.DOC .* g.V)+sum(gtr.POC .* g.V)
    ΣsurFⁿ= sum((nutₜ.DIN[:,:,1]+nutₜ.DON[:,:,1]+nutₜ.PON[:,:,1]) .* g.Az .* vel.w[:,:,1])
    ΣsurFᶜ= sum((nutₜ.DIC[:,:,1]+nutₜ.DOC[:,:,1]+nutₜ.POC[:,:,1]) .* g.Az .* vel.w[:,:,1])
    ΣDIN = sum(nutₜ.DIN .* g.V)
    Cio = open("results/cons_C.txt","a"); Nio = open("results/cons_N.txt","a");
    DINio = open("results/cons_DIN.txt","a");
    println(Cio,@sprintf("%3.0f  %.16E  %.16E  %.8E",t,Σgtrᶜ,ΣsurFᶜ,Σgtrᶜ+ΣsurFᶜ))
    println(Nio,@sprintf("%3.0f  %.16E  %.16E  %.8E",t,Σgtrⁿ,ΣsurFⁿ,Σgtrⁿ+ΣsurFⁿ))
    println(DINio,@sprintf("%3.0f  %.16E %7.0f %5.0f %5.0f %5.0f",t,ΣDIN,agent_num,death_ct,graz_ct,dvid_ct))
    close(Cio);close(Nio);close(DINio);
end

"""
    PrepRunDir(res::String="results/")
Create `res/` folder if needed. Remove old files from it if needed.
"""
function PrepRunDir(res::String="results/")
    isdir(res) || mkdir(res)
    isdir("$res"*"nutrients/") || mkdir("$res"*"nutrients/")
    println("$res"*"nutrients/")

    isfile("$res"*"cons_C.txt") && rm("$res"*"cons_C.txt");
    isfile("$res"*"cons_N.txt") && rm("$res"*"cons_N.txt");
    isfile("$res"*"cons_DIN.txt") && rm("$res"*"cons_DIN.txt");
    isfile("$res"*"B1.bin") && rm("$res"*"B1.bin");
    isfile("$res"*"B2.bin") && rm("$res"*"B2.bin");
    isfile("$res"*"output.bin") && rm("$res"*"output.bin");
    isfile("$res"*"output1.bin") && rm("$res"*"output1.bin");
    isfile("$res"*"output2.bin") && rm("$res"*"output2.bin");
    isfile("$res"*"grid.bin") && rm("$res"*"grid.bin");
    isfile("$res"*"IR.bin") && rm("$res"*"IR.bin");
    isfile("$res"*"VD1.bin") && rm("$res"*"VD1.bin");
    isfile("$res"*"VD2.bin") && rm("$res"*"VD2.bin");
    isfile("$res"*"HD1.bin") && rm("$res"*"HD1.bin");
    isfile("$res"*"HD2.bin") && rm("$res"*"HD2.bin");

    return "done"

end

"""
  testB1B2(B1,B2,fil="")
Compute mean size for `B1` & `B2`, add `B1ref` & `B2ref` reference result obtained from `fil`, and output all time series in one DataFrame.

```
tst=testB1B2(B1,B2,"samples/testB1B2.csv")
tst2=[tst[!,:B1],tst[!,:B1ref],tst[!,:B2],tst[!,:B2ref]]
isapprox(tst[end,:B1],tst[end,:B1ref]; atol=1e-2)
isapprox(tst[end,:B2],tst[end,:B2ref]; atol=1e-2)
using Plots; plot(tst2,lab=["B1" "B1ref" "B2" "B2ref"])
```
"""
function testB1B2(B1,B2,fil="")
    n=length(B1)
    tmpB1=map(x -> mean(B1[x][!,:size]), 1:n)
    tmpB2=map(x -> mean(B2[x][!,:size]), 1:n)
    df=DataFrame(B1=tmpB1,B2=tmpB2)
    if ~isempty(fil)
        tmpB1B2=CSV.read(fil)
        df[!,:B1ref]=tmpB1B2[!,:B1]
        df[!,:B2ref]=tmpB1B2[!,:B2]
    end
    return df
end
