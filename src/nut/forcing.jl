@kernel function calc_forcing!(F, nutrients, params, ΔT)
    i, j, k = @index(Global, NTuple)
    F.DIC.data[i, j, k] = max(0.0, nutrients.DOC.data[i, j, k]) * params["kDOC"] * ΔT
    F.DOC.data[i, j, k] = max(0.0, nutrients.POC.data[i, j, k]) * params["kPOC"] * ΔT -
                          max(0.0, nutrients.DOC.data[i, j, k]) * params["kDOC"] * ΔT
    F.POC.data[i, j, k] =-max(0.0, nutrients.POC.data[i, j, k]) * params["kPOC"] * ΔT

    F.NH4.data[i, j, k] = max(0.0, nutrients.DON.data[i, j, k]) * params["kDON"] * ΔT -
                          max(0.0, nutrients.NH4.data[i, j, k]) * params["Nit"]  * ΔT
    F.NO3.data[i, j, k] = max(0.0, nutrients.NH4.data[i, j, k]) * params["Nit"]  * ΔT
    F.DON.data[i, j, k] = max(0.0, nutrients.PON.data[i, j, k]) * params["kPON"] * ΔT -
                          max(0.0, nutrients.DON.data[i, j, k]) * params["kDON"] * ΔT
    F.PON.data[i, j, k] =-max(0.0, nutrients.PON.data[i, j, k]) * params["kPON"] * ΔT

    F.PO4.data[i, j, k] = max(0.0, nutrients.DOP.data[i, j, k]) * params["kDOP"] * ΔT
    F.DOP.data[i, j, k] = max(0.0, nutrients.POP.data[i, j, k]) * params["kPOP"] * ΔT -
                          max(0.0, nutrients.DOP.data[i, j, k]) * params["kDOP"] * ΔT
    F.POP.data[i, j, k] =-max(0.0, nutrients.POP.data[i, j, k]) * params["kPOP"] * ΔT
end

function nut_forcing!(F, arch::Architecture, g, nutrients, params, ΔT)
    calc_forcing_kernel! = calc_forcing!(device(arch), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))

    event = calc_forcing_kernel!(F, nutrients, params, ΔT)

    wait(device(arch), event)

    return nothing
end

"""
    sub_nut_tendency!(a, b, c)
subtract one tendency from total tendencies
"""
function sub_nut_tendency!(a, b, c)
    for i in 1:length(a)
        a[i].data .= b[i].data .- c[i].data
    end
end

"""
    add_nut_tendency!(a, b)
add one tendency to another tendency
"""
function add_nut_tendency!(a, b)
    for i in 1:length(a)
        a[i].data .= a[i].data .+ b[i].data
    end
end
