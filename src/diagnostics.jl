#=
Diagnostics of individual physiology processes
Indices:
1:  Counts, number of individuals appear in each grid during the time between two diagnostics
2:  PP
3:  VNO3
4:  VHN4
5:  VP
6:  VDOC
7:  BS
8:  Respiration
9:  Exudation
10: Bm
11: Cq
12: Nq
13: Pq
14: Chl
=#
"""
"""
function diags_setup(nTime::Int64, ΔT::Int64, grids, freq::Int64, diag_inds::Array, Nsp::Int64)
    ndiags = sum(diag_inds)
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags_sp = zeros(grids.Nx, grids.Ny, grids.Nz, nt, Nsp, ndiags)
    diags_pop = zeros(grids.Nx, grids.Ny, grids.Nz, nt, Nsp, 3)
    return diags_sp, diags_pop
end

function diags_setup(nTime::Int64, ΔT::Int64, grids, freq::Int64, nTr::Int64)
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags = zeros(grids.Nx, grids.Ny, grids.Nz, nt, nTr)
    return diags
end
