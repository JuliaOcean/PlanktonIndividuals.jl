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
