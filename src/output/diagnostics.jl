#=
Diagnostics of individual physiology processes
Indices:
1:  counts
2:  grazing
3:  mortality
4:  division
5:  PS
6:  VDOC
7:  VNO3
8:  VHN4
9:  VPO4
10: Respiration
11: BS
12: Exudation
13: Bm
14: Cq
15: Nq
16: Pq
17: Chl
=#
mutable struct Diagnostics
    spcs::AbstractArray{Float64,6}       # for each species
    tr::AbstractArray{Float64,5}         # for tracers
end

function diags_setup(arch::Architecture, nTime::Int64, ΔT::Int64, grids, freq::Int64, Nsp::Int64, nTr::Int64)
    ndiags = 17
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags_sp = zeros(grids.Nx, grids.Ny, grids.Nz, nt, Nsp, ndiags) |> array_type(arch)
    diags = zeros(grids.Nx, grids.Ny, grids.Nz, nt, nTr) |> array_type(arch)
    return Diagnostics(diags_sp, diags)
end


##### record diagnostics at each time step
@kernel function sum_diags_kernel!(diags, plank, inds::AbstractArray{Int64,2},
                                   sp::Int64, g::Grids, diag_t)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]

        @inbounds diags[xi, yi, zi, diag_t, sp, 1]  += 1 # individual count
        @inbounds diags[xi, yi, zi, diag_t, sp, 2]  += plank[i,31] # grazing

        @inbounds diags[xi, yi, zi, diag_t, sp, 5]  += plank[i,22] # PS
        @inbounds diags[xi, yi, zi, diag_t, sp, 6]  += plank[i,23] # VDOC
        @inbounds diags[xi, yi, zi, diag_t, sp, 7]  += plank[i,24] # VNH4
        @inbounds diags[xi, yi, zi, diag_t, sp, 8]  += plank[i,25] # VNO3
        @inbounds diags[xi, yi, zi, diag_t, sp, 9]  += plank[i,26] # VPO4
        @inbounds diags[xi, yi, zi, diag_t, sp, 10] += plank[i,28] # respir
        @inbounds diags[xi, yi, zi, diag_t, sp, 11] += plank[i,29] # BS
        @inbounds diags[xi, yi, zi, diag_t, sp, 12] += plank[i,30] # exu
        @inbounds diags[xi, yi, zi, diag_t, sp, 13] += plank[i,6]  # Bm
        @inbounds diags[xi, yi, zi, diag_t, sp, 14] += plank[i,7]  # Cq
        @inbounds diags[xi, yi, zi, diag_t, sp, 15] += plank[i,8]  # Nq
        @inbounds diags[xi, yi, zi, diag_t, sp, 16] += plank[i,9]  # Pq
        @inbounds diags[xi, yi, zi, diag_t, sp, 17] += plank[i,10] # Chl
    end
end
function sum_diags!(diags, plank, inds::AbstractArray{Int64,2},
                    sp::Int64, arch::Architecture, g::Grids, diag_t)
    kernel! = sum_diags_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(diags, plank, inds, sp, g, diag_t)
    wait(device(arch), event)
    return nothing
end

@kernel function sum_diags_mort_kernel!(diags, plank, inds::AbstractArray{Int64,2},
                                        sp::Int64, g::Grids, diag_t)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds diags[xi, yi, zi, diag_t, sp, 3] += plank[i,32]
    end
end
function sum_diags_mort!(diags, plank, inds::AbstractArray{Int64,2},
                         sp::Int64, arch::Architecture, g::Grids, diag_t)
    kernel! = sum_diags_mort_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(diags, plank, inds, sp, g, diag_t)
    wait(device(arch), event)
    return nothing
end

@kernel function sum_diags_dvid_kernel!(diags, plank, inds::AbstractArray{Int64,2},
                                        sp::Int64, g::Grids, diag_t)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds diags[xi, yi, zi, diag_t, sp, 4] += plank[i,33]
    end
end
function sum_diags_dvid!(diags, plank, inds::AbstractArray{Int64,2},
                         sp::Int64, arch::Architecture, g::Grids, diag_t)
    kernel! = sum_diags_dvid_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(diags, plank, inds, sp, g, diag_t)
    wait(device(arch), event)
    return nothing
end

