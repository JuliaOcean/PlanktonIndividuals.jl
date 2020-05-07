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
8:  Mainten
9:  Exudation
10: Division
11: Grazing
12: Death
13: Bm
14: Cq
15: Nq
16: Pq
17: Chl
18: PAR
=#
"""
"""
function diags_setup(nTime::Int64, ΔT::Int64, grids, freq::Int64, diag_inds::Array)
    ndiags = sum(diag_inds)
    nt = nTime*ΔT÷freq
    diags = zeros(grids.Nx, grids.Ny, grids.Nz,nt,ndiags)
    return diags
end
