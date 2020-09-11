@kernel function calc_DIC_forcing!(F_DIC, grid, DOC, kDOC, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_DIC[ii, jj, kk] = F_DIC[ii, jj, kk] + max(0.0, DOC[ii, jj, kk]) * kDOC * ΔT
end

@kernel function calc_DOC_forcing!(F_DOC, grid, DOC, POC, kDOC, kPOC, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_DOC[ii, jj, kk] = F_DOC[ii, jj, kk] + max(0.0, POC[ii, jj, kk]) * kPOC * ΔT -
                                                      max(0.0, DOC[ii, jj, kk]) * kDOC * ΔT
end

@kernel function calc_POC_forcing!(F_POC, grid, POC, kPOC, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_POC[ii, jj, kk] =F_POC[ii, jj, kk] - max(0.0, POC[ii, jj, kk]) * kPOC * ΔT
end

@kernel function calc_NH4_forcing!(F_NH4, grid, NH4, DON, kNit, kDON, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_NH4[ii, jj, kk] = F_NH4[ii, jj, kk] + max(0.0, DON[ii, jj, kk]) * kDON * ΔT -
                                                      max(0.0, NH4[ii, jj, kk]) * kNit * ΔT
end

@kernel function calc_NO3_forcing!(F_NO3, grid, NH4, kNit, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_NO3[ii, jj, kk] =  F_NO3[ii, jj, kk] + max(0.0, NH4[ii, jj, kk]) * kNit * ΔT
end

@kernel function calc_DON_forcing!(F_DON, grid, DON, PON, kDON, kPON, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_DON[ii, jj, kk] = F_DON[ii, jj, kk] + max(0.0, PON[ii, jj, kk]) * kPON * ΔT -
                                                      max(0.0, DON[ii, jj, kk]) * kDON * ΔT
end

@kernel function calc_PON_forcing!(F_PON, grid, PON, kPON, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_PON[ii, jj, kk] = F_PON[ii, jj, kk] - max(0.0, PON[ii, jj, kk]) * kPON * ΔT
end

@kernel function calc_PO4_forcing!(F_PO4, grid, DOP, kDOP, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_PO4[ii, jj, kk] = F_PO4[ii, jj, kk] + max(0.0, DOP[ii, jj, kk]) * kDOP * ΔT
end

@kernel function calc_DOP_forcing!(F_DOP, grid, DOP, POP, kDOP, kPOP, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_DOP[ii, jj, kk] = F_DOP[ii, jj, kk] + max(0.0, POP[ii, jj, kk]) * kPOP * ΔT -
                                                      max(0.0, DOP[ii, jj, kk]) * kDOP * ΔT
end

@kernel function calc_POP_forcing!(F_POP, grid, POP, kPOP, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds F_POP[ii, jj, kk] = F_POP[ii, jj, kk] - max(0.0, POP[ii, jj, kk]) * kPOP * ΔT
end

function nut_forcing!(F, arch::Architecture, g, nut, params, ΔT)
    calc_DIC_forcing_kernel! = calc_DIC_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_DOC_forcing_kernel! = calc_DOC_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_POC_forcing_kernel! = calc_POC_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_NH4_forcing_kernel! = calc_NH4_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_NO3_forcing_kernel! = calc_NO3_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_DON_forcing_kernel! = calc_DON_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_PON_forcing_kernel! = calc_PON_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_PO4_forcing_kernel! = calc_PO4_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_DOP_forcing_kernel! = calc_DOP_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    calc_POP_forcing_kernel! = calc_POP_forcing!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))

    ev_DIC = calc_DIC_forcing_kernel!(F.DIC.data, g, nut.DOC.data, params["kDOC"], ΔT)
    ev_DOC = calc_DOC_forcing_kernel!(F.DOC.data, g, nut.DOC.data, nut.POC.data, params["kDOC"], params["kPOC"], ΔT)
    ev_POC = calc_POC_forcing_kernel!(F.POC.data, g, nut.POC.data, params["kPOC"], ΔT)
    ev_NH4 = calc_NH4_forcing_kernel!(F.NH4.data, g, nut.NH4.data, nut.DON.data, params["Nit"], params["kDON"], ΔT)
    ev_NO3 = calc_NO3_forcing_kernel!(F.NO3.data, g, nut.NH4.data, params["Nit"], ΔT)
    ev_DON = calc_DON_forcing_kernel!(F.DON.data, g, nut.DON.data, nut.PON.data, params["kDON"], params["kPON"], ΔT)
    ev_PON = calc_PON_forcing_kernel!(F.PON.data, g, nut.PON.data, params["kPON"], ΔT)
    ev_PO4 = calc_PO4_forcing_kernel!(F.PO4.data, g, nut.DOP.data, params["kDOP"], ΔT)
    ev_DOP = calc_DOP_forcing_kernel!(F.DOP.data, g, nut.DOP.data, nut.POP.data, params["kDOP"], params["kPOP"], ΔT)
    ev_POP = calc_POP_forcing_kernel!(F.POP.data, g, nut.POP.data, params["kPOP"], ΔT)

    events = [ev_DIC, ev_DOC, ev_POC, ev_NH4, ev_NO3, ev_DON, ev_PON, ev_PO4, ev_DOP, ev_POP]

    wait(device(arch), MultiEvent(Tuple(events)))

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

function zero_fields!(a)
    for i in 1:length(a)
        a[i].data .= 0.0
    end
end
# """
#     add_nut_tendency!(a, b)
# add one tendency to another tendency
# """
# function add_nut_tendency!(a, b)
#     for i in 1:length(a)
#         a[i].data .= a[i].data .+ b[i].data
#     end
# end
