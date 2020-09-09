#=
Diagnostics of individual physiology processes
Indices:
1:  PP
2:  VDOC
3:  VNO3
4:  VHN4
5:  VP
6:  Respiration
7:  BS
8:  Exudation
9:  Bm
10: Cq
11: Nq
12: Pq
13: Chl
=#
mutable struct Diagnostics
    spcs::AbstractArray{Float64,6}       # for each species
    tr::AbstractArray{Float64,5}         # for tracers
end

function diags_setup(::CPUs, nTime::Int64, ΔT::Int64, grids, freq::Int64, diag_inds::Array, Nsp::Int64)
    ndiags = sum(diag_inds) + 4
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags_sp = zeros(grids.Nx, grids.Ny, grids.Nz, nt, Nsp, ndiags)
    return diags_sp
end
function diags_setup(::GPUs, nTime::Int64, ΔT::Int64, grids, freq::Int64, diag_inds::Array, Nsp::Int64)
    ndiags = sum(diag_inds) + 4
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags_sp = zeros(grids.Nx, grids.Ny, grids.Nz, nt, Nsp, ndiags) |> CuArray
    return diags_sp
end

function diags_setup(::CPUs, nTime::Int64, ΔT::Int64, grids, freq::Int64, nTr::Int64)
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags = zeros(grids.Nx, grids.Ny, grids.Nz, nt, nTr)
    return diags
end
function diags_setup(::GPUs, nTime::Int64, ΔT::Int64, grids, freq::Int64, nTr::Int64)
    nt = nTime*ΔT÷freq
    if nTime*ΔT%freq ≠ 0
        nt += 1
    end
    diags = zeros(grids.Nx, grids.Ny, grids.Nz, nt, nTr) |> CuArray
    return diags
end

##### record diagnostics at each time step
@kernel function sum_diags_kernel!(diags, op_array, g::Grids, diags_inds, diag_t)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[i,1], g) |> Int
    yi = find_yF_ind(op_array[i,2], g) |> Int
    zi = find_zF_ind(op_array[i,3], g) |> Int
    sp = op_array[i,11]

    diags[xi, yi, zi, diag_t, sp, 1] += 1 # individual count
    diags[xi, yi, zi, diag_t, sp, 2] += op_array[i,29] # grazing

    id = 4
    if diags_inds[1] == 1 # PS
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,20]
    end
    if diags_inds[2] == 1 # VDOC
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,21]
    end
    if diags_inds[3] == 1 # VNH4
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,22]
    end
    if diags_inds[4] == 1 # VNO3
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,23]
    end
    if diags_inds[5] == 1 # VPO4
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,24]
    end
    if diags_inds[6] == 1 # respiration
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,26]
    end
    if diags_inds[7] == 1 # BS
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,27]
    end
    if diags_inds[8] == 1 # exudation
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,28]
    end
    if diags_inds[9] == 1 # Bm
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,6]
    end
    if diags_inds[10] == 1 # Cq
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,7]
    end
    if diags_inds[11] == 1 # Nq
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,8]
    end
    if diags_inds[12] == 1 # Pq
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,9]
    end
    if diags_inds[13] == 1 # Chla
        id += 1
        diags[xi, yi, zi, diag_t, sp, id] += op_array[i,10]
    end
end
function sum_diags!(diags, op_array, arch::Architecture, g::Grids, diags_inds, diag_t)
    kernel! = sum_diags_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(diags, op_array, g, diags_inds, diag_t)
    wait(device(arch), event)
    return nothing
end

@kernel function sum_diags_mort_kernel!(diags_spcs, op_array, g::Grids, diag_t)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[i,1], g) |> Int
    yi = find_yF_ind(op_array[i,2], g) |> Int
    zi = find_zF_ind(op_array[i,3], g) |> Int
    sp = op_array[i,11]
    diags[xi, yi, zi, diag_t, sp, 3] += op_array[i,30]
end
function sum_diags_mort!(diags, op_array, arch::Architecture, g::Grids, diag_t)
    kernel! = sum_diags_mort_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(diags, op_array, g, diag_t)
    wait(device(arch), event)
    return nothing
end

@kernel function sum_diags_dvid_kernel!(diags_spcs, op_array, g::Grids, diag_t)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[i,1], g) |> Int
    yi = find_yF_ind(op_array[i,2], g) |> Int
    zi = find_zF_ind(op_array[i,3], g) |> Int
    sp = op_array[i,11]
    diags[xi, yi, zi, diag_t, sp, 4] += op_array[i,31]
end
function sum_diags_dvid!(diags, op_array, arch::Architecture, g::Grids, diag_t)
    kernel! = sum_diags_dvid_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(diags, op_array, g, diag_t)
    wait(device(arch), event)
    return nothing
end


