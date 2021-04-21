"""
    update_bgc_params(tmp::Ditc)
Update parameter values based on a `Dict` provided by user

Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
"""
function update_bgc_params(tmp::Dict)
    parameters = bgc_params_default()
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: parameter not found $key"))
        else
            parameters[key] = tmp[key]
        end
    end
    return parameters
end

"""
    update_phyt_params(tmp::Dict)
Update parameter values based on a `Dict` provided by user

Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
"""
function update_phyt_params(tmp::Dict)
    parameters = phyt_params_default()
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: parameter not found $key"))
        else
            parameters[key] = tmp[key]
        end
    end
    return parameters
end

function read_IR_input(ΔT::Int64, grid; path)
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

function read_temp_input(ΔT::Int64, grid; path, ∂T∂z=0.04)
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
    return temp_domain
end

