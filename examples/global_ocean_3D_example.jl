
import MeshArrays, OceanStateEstimation, NetCDF

MeshArrays.GRID_LL360_download()
γ = MeshArrays.GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
Γ = MeshArrays.GridLoad(γ;option="full")

grid_info=( XC=Γ.XC[1], YC=Γ.YC[1], RC=Γ.RC,
            XW=Γ.XW[1], YS=Γ.YS[1], RF=Γ.RF,
            DXC=Γ.DXC[1], DYC=Γ.DYC[1], DRC=Γ.DRC,
            DXG=Γ.DXG[1], DYG=Γ.DYG[1], DRF=Γ.DRF,
            RAW=Γ.RAW[1], RAS=Γ.RAS[1], RAC=Γ.RAC[1],
            hFacC=write(Γ.hFacC), hFacW=write(Γ.hFacW), hFacS=write(Γ.hFacS));

OceanStateEstimation.get_occa_variable_if_needed("DDtheta")
OceanStateEstimation.get_occa_variable_if_needed("DDuvel")
OceanStateEstimation.get_occa_variable_if_needed("DDvvel")
OceanStateEstimation.get_occa_variable_if_needed("DDwvel")

function rd(filename, varname, month)
    fil = NetCDF.open(filename, varname)
    tmp = fil[:,:,:,month]
    tmp[findall(tmp.<-1e22)] .= 0.0
    return tmp
end

month=1
pth=OceanStateEstimation.OCCAclim_path
U=rd(joinpath(pth,"DDuvel.0406clim.nc"),"u",month);
V=rd(joinpath(pth,"DDvvel.0406clim.nc"),"v",month);
W=rd(joinpath(pth,"DDwvel.0406clim.nc"),"w",month);
T=rd(joinpath(pth,"DDtheta.0406clim.nc"),"theta",month);
