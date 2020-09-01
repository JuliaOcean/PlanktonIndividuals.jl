@kernel function calc_forcing!(F, nutrients, params, ΔT)
    i, j, k = @index(Global, NTuple)
    F.DIC[i, j, k] = max(0.0, nutrients.DOC[i, j, k]) * params["kDOC"] * ΔT
    F.DOC[i, j, k] = max(0.0, nutrients.POC[i, j, k]) * params["kPOC"] * ΔT -
                     max(0.0, nutrients.DOC[i, j, k]) * params["kDOC"] * ΔT
    F.POC[i, j, k] =-max(0.0, nutrients.POC[i, j, k]) * params["kPOC"] * ΔT

    F.NH4[i, j, k] = max(0.0, nutrients.DON[i, j, k]) * params["kDON"] * ΔT -
                     max(0.0, nutrients.NH4[i, j, k]) * params["Nit"]  * ΔT
    F.NO3[i, j, k] = max(0.0, nutrients.NH4[i, j, k]) * params["Nit"]  * ΔT
    F.DON[i, j, k] = max(0.0, nutrients.PON[i, j, k]) * params["kPON"] * ΔT -
                     max(0.0, nutrients.DON[i, j, k]) * params["kDON"] * ΔT
    F.PON[i, j, k] =-max(0.0, nutrients.PON[i, j, k]) * params["kPON"] * ΔT

    F.PO4[i, j, k] = max(0.0, nutrients.DOP[i, j, k]) * params["kDOP"] * ΔT
    F.DOP[i, j, k] = max(0.0, nutrients.POP[i, j, k]) * params["kPOP"] * ΔT -
                     max(0.0, nutrients.DOP[i, j, k]) * params["kDOP"] * ΔT
    F.POP[i, j, k] =-max(0.0, nutrients.POP[i, j, k]) * params["kPOP"] * ΔT
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
function sub_nut_tendency!(a::nutrient_fields, b::nutrient_fields, c::nutrient_fields)
    a.DOC = b.DOC .- c.DOC
    a.POC = b.POC .- c.POC
    a.DON = b.DON .- c.DON
    a.PON = b.PON .- c.PON
    a.DOP = b.DOP .- c.DOP
    a.POP = b.POP .- c.POP
    a.DIC = b.DIC .- c.DIC
    a.NH4 = b.NH4 .- c.NH4
    a.NO3 = b.NO3 .- c.NO3
    a.PO4 = b.PO4 .- c.PO4
end
"""
    add_nut_tendency!(a, b)
add one tendency to another tendency
"""
function add_nut_tendency!(a::nutrient_fields, b::nutrient_fields)
    a.DOC = a.DOC .+ b.DOC
    a.POC = a.POC .+ b.POC
    a.DON = a.DON .+ b.DON
    a.PON = a.PON .+ b.PON
    a.DOP = a.DOP .+ b.DOP
    a.POP = a.POP .+ b.POP
    a.DIC = a.DIC .+ b.DIC
    a.NH4 = a.NH4 .+ b.NH4
    a.NO3 = a.NO3 .+ b.NO3
    a.PO4 = a.PO4 .+ b.PO4
end
