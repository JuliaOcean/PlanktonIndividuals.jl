"""
    update_params!(parameters, tmp)
Update parameter values based on a .yaml file provided by user

Keyword Arguments
=================
- `parameters` is default parameter set
- `tmp` is the parameters read from .yaml file and needs to update
"""
function update_params!(parameters::Dict, tmp::Dict)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: parameter not found $key"))
        else
            parameters[key] = tmp[key]
        end
    end
end

function read_IR_input(ΔT::Int64,grid,
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
    PAR = itp_PAR.(t_ΔT)
    # expand to the whole domain surface
    PAR_domain = zeros(grid.Nx, grid.Ny, size(PAR,1))
    for i in 1:size(PAR,1)
        PAR_domain[:,:,i] .= PAR[i]
    end
    return PAR_domain
end

function read_temp_input(ΔT::Int64, grid, ∂T∂z=0.04,
                         path = dirname(pathof(PlanktonIndividuals))*"/../samples/temp.bin")
    temp_hour = deserialize(path)
    # convert hour to second in a day
    t_htos = collect(0:1:24) .* 3600
    # convert each time step to second in a day
    t_ΔT = collect(1:ΔT:86400)
    # interpolation of PAR time series
    itp_temp = interpolate((t_htos,), temp_hour, Gridded(Linear()));
    temp = itp_temp.(t_ΔT)
    # espand to the whole domain
    temp_domain = zeros(grid.Nx, grid.Ny, grid.Nz, size(temp,1))
    for i in 1:size(temp,1)
        temp_domain[:,:,end,i] .= temp[i]
    end
    # vertical temperature gradient
    for j in grid.Nz-1:-1:1
        temp_domain[:,:,j,:] .= temp_domain[:,:,j+1,:] .- (∂T∂z*(grid.zC[j+1]-grid.zC[j]))
    end
    temp_d = zeros(grid.Nx+grid.Hx*2, grid.Ny+grid.Hy*2, grid.Nz+grid.Hz*2, size(temp,1))
    temp_d[grid.Hx+1:grid.Nx+grid.Hx, grid.Hy+1:grid.Ny+grid.Hy, grid.Hz+1:grid.Nz+grid.Hz, :] .= temp_domain
    return temp_d
end

