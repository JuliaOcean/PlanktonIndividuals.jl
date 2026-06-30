"""
    update_bgc_params(tmp::Ditc, FT::DataType)
Update parameter values based on a `Dict` provided by user

Keyword Arguments
=================
- `tmp`: a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
"""
function update_bgc_params(tmp::Dict, FT::DataType)
    parameters = bgc_params_default(FT)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: bgc parameter not found $key"))
        else
            parameters[key] = FT(tmp[key])
        end
    end
    return parameters
end

"""
    update_phyt_params(tmp::Dict, FT::DataType; N::Int64, mode::AbstractMode)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
- `mode` is the mode of phytoplankton physiology resolved in the model
"""
function update_phyt_params(tmp::Dict, FT::DataType; N::Int = 1, mode::AbstractMode = QuotaMode())
    parameters = phyt_params_default(N,mode)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: phyt parameter not found $key"))
        else
            parameters[key] = FT.(tmp[key])
        end
    end
    return parameters
end

"""
    update_colony_params(tmp::Dict, FT::DataType; N::Int64, mode::AbstractMode)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
- `mode` is the mode of phytoplankton physiology resolved in the model
"""
function update_colony_params(tmps::AbstractArray, FT::DataType; 
                              Ncl::Int = 1, Nsp::AbstractArray = [1], 
                              mode::AbstractMode = IronEnergyMode())
    parameters = colony_params_default(Ncl, Nsp, mode)
    for i in eachindex(tmps)
        tmp = tmps[i]
        tmp_keys = collect(keys(tmp))
        pkeys = collect(keys(parameters[i]))
        for key in tmp_keys
            if length(findall(x->x==key, pkeys))==0
                throw(ArgumentError("PARAM: colony parameter not found $key"))
            else
                parameters[i][key] = FT.(tmp[key])
            end
        end
    end
    return parameters
end


"""
    update_abiotic_params(tmp::Dict, FT::DataType; N::Int64)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
"""
function update_abiotic_params(tmp::Dict, FT::DataType; N::Int = 1)
    parameters = abiotic_params_default(N)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: abiotic parameter not found $key"))
        else
            parameters[key] = FT.(tmp[key])
        end
    end
    return parameters
end
