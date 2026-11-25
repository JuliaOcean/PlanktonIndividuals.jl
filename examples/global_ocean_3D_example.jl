# # Global Ocean 3D example
#
# Here we setup a passive tracer and particle simulation. The tracer is initialized to 1 below 100m depth and 0 above that depth horizon to highlight the Eulerian advection capability. 1024 particles are initially distributed at random over the Global Ocean to demonstrate the Lagrangian advection capability.
#
# !!! note
#     This simulation may take several hours to complete a one-year run. Here we only run a few time steps.

using PlanktonIndividuals, PlanktonKernels
import Downloads, MeshArrays, NetCDF, Plots, Scratch

# Load the grid and cache the [OCCA climatology](https://doi.org/10.7910/DVN/RNXA2A) locally.

begin
    MeshArrays.GRID_LL360_download();
    γ = MeshArrays.GridSpec("PeriodicChannel",MeshArrays.GRID_LL360);
    Γ = MeshArrays.GridLoad(γ;option="full");

    grid_info=(XC=Γ.XC[1], YC=Γ.YC[1], RC=Γ.RC,
            XW=Γ.XW[1], YS=Γ.YS[1], RF=Γ.RF,
            DXC=Γ.DXC[1], DYC=Γ.DYC[1], DRC=Γ.DRC,
            DXG=Γ.DXG[1], DYG=Γ.DYG[1], DRF=Γ.DRF,
            RAW=Γ.RAW[1], RAS=Γ.RAS[1], RAC=Γ.RAC[1],
            hFacC=ones(360,160,50), hFacW=ones(360,160,50), hFacS=ones(360,160,50));

    pth = Scratch.get_scratch!(PlanktonIndividuals, "occa_climatology")
    for (name, id) in (("DDtheta", 2842852), ("DDuvel", 2842849),
                       ("DDvvel", 2842857), ("DDwvel", 2842856))
        file = joinpath(pth, name * ".0406clim.nc")
        if !isfile(file)
            downloaded = Downloads.download("https://dataverse.harvard.edu/api/access/datafile/$id?format=original")
            mv(downloaded, file; force=true)
        end
    end
    nothing
end

function rd(filename, varname, month)
    values = NetCDF.ncread(filename, varname; start=[1, 1, 1, month], count=[-1, -1, -1, 1])
    tmp = dropdims(values; dims=4)
    tmp[tmp .< -1e22] .= 0.0
    return tmp
end

# The five hourly steps use only the first monthly field.

begin
    u = Float32.(rd(joinpath(pth, "DDuvel.0406clim.nc"), "u", 1))
    v = zeros(Float32, 360, 161, 50, 1)
    w = zeros(Float32, 360, 160, 51, 1)
    v[:, 1:160, :, 1] .= rd(joinpath(pth, "DDvvel.0406clim.nc"), "v", 1)
    w[:, :, 2:51, 1] .= rd(joinpath(pth, "DDwvel.0406clim.nc"), "w", 1)
    temperature = Float32.(rd(joinpath(pth, "DDtheta.0406clim.nc"), "theta", 1))
    uvels = reshape(u, size(u)..., 1)
    vvels, wvels = v, w
    temps = reshape(temperature, size(temperature)..., 1)
    mask = Array(temperature .!= 0)
    nothing
end

grid = LoadLatLonGrid(;grid_info = grid_info, size = (360,160,50), lat = (-80,80), lon= (-180,180), landmask = mask)

phyt_parameter = Dict("PCmax"   => [0.0], # maximum photosynthesis rate
                      "respir"=> [0.0], # respiration rate
                      "dvid_P"  => [0.0], # probability of cell division
                      "mort_P"  => [0.0], # probability of cell natural death
                      "thermal" => [0.0]  # thermal damage and repair
                     )

bgc_parameter = Dict("kDOC" => 0.0, # Remineralization rate for DOC
                     "Nit"  => 0.0, # Nitrification rate for NH4
                     "kDON" => 0.0, # Remineralization rate for DON
                     "kDOP" => 0.0, # Remineralization rate for DON
                     "kPOC" => 0.0, # Remineralization rate for POC
                     "kPON" => 0.0, # Remineralization rate for PON
                     "kPOP" => 0.0, # Remineralization rate for POP
                     "w_sink_org" => 0.0,
                     "w_sink_inorg" => 0.0,
                    )

model = PlanktonModel(CPU(), grid;
                      mode = CarbonMode(),
                      phyto = phyto_setup(phyt_parameter, [1024], 1),
                      options = model_options(max_individuals = 2^13),
                      bgc_params = bgc_parameter)

begin
    tmp=1.0*mask
    tmp[:,:,1:10].=0.0
    DOC_init=zeros(size(mask).+4)
    DOC_init[3:end-2,3:end-2,3:end-2] .= tmp
    model.tracers.DOC.data.=DOC_init
    ✔2= "Done with initial condition"
end

begin
    ✔2
    sim = PlanktonSimulation(model, ΔT = 3600.0,
                                    iterations = 5,
                                    vels=(u=uvels, v=vvels, w=wvels),
                                    temp = temps,
                                    ΔT_vel=3600.0*24*30,
                                    ΔT_temp=3600.0*24*30)
end

simulated_model = begin
    sim_run!(sim)
    sim.model
end

function plot_map(model::PlanktonModel, kk=10)
    ## tracer field
    col=Array(model.tracers.DOC.data)[3:end-2,3:end-2,kk]
    col[findall(mask[:,:,kk].==0.0)].=NaN
    xC, yC = collect(model.grid.xC)[3:end-2], collect(model.grid.yC)[3:end-2]
    fig = Plots.heatmap(xC, yC, col',clims=(-0.1, 1.1), color=:thermal, fmt=:png)

    ## scatter plot
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 .- 180 # convert fractional indices to degree
    py = Array(model.individuals.phytos.sp1.data.y) .* 1 .- 80  # convert fractional indices to degree
    Plots.scatter!(fig, px, py, ms=3, color = :white, legend=:none)

    return fig
end

# The map shows the simulated tracer and particle positions.

f1 = plot_map(simulated_model, 3)

# The zonal mean shows how the tracer varies with latitude and depth.

begin
    function plot_zm(model::PlanktonModel, kk=10)
        d1=Array(model.tracers.DOC.data)[3:end-2,3:end-2,3:end-2];
        #d1[findall(mask.==0.0)].=missing
        d2=dropdims(sum(d1,dims=1),dims=1)./dropdims(sum(mask,dims=1),dims=1)
        zC2, yC2 = collect(model.grid.zC)[3:end-2], collect(model.grid.yC)[3:end-2]
        f2=Plots.contourf(yC2,reverse(zC2),reverse(permutedims(d2),dims=1),
        ylims=(-200.0,0.0),xlabel="latitude",ylabel="depth")
        f2
    end

    f2=plot_zm(simulated_model)
end
