mutable struct Diagnostics
    spcs::NamedTuple       # for each species
    tr::NamedTuple         # for tracers
end

function diags_setup(ntrs, nprocs, g::RegularRectilinearGrid, Nsp::Int64, arch::Architecture)
    ntr   = length(ntrs)
    nproc = length(nprocs)
    trs   = []
    procs = []

    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)

    for i in 1:ntr
        tr = zeros(total_size) |> array_type(arch)
        push!(trs, tr)
    end

    diag_tr = NamedTuple{ntrs}(trs)

    for j in 1:Nsp
        procs_sp = []
        for k in 1:nproc
            proc = zeros(total_size) |> array_type(arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{nprocs}(procs_sp)
        push!(procs, diag_proc)
    end
    plank_name = plank_names[1:Nsp]
    diag_sp = NamedTuple{plank_name}(procs)
    return Diagnostics(diag_sp, diag_tr)
end

##### record diagnostics at each time step
function gpu_diags_proc_kernel!(diags_proc, proc, ac, x, y, z)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds @atomic diags_proc[x[i], y[i], z[i]] += proc[i] * ac[i]
    end
end
function diags_proc!(diags_proc, proc, ac, x, y, z, ::GPUs)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_diags_proc_kernel!(diags_proc, proc, ac, x, y, z)
    return nothing 
end
function diags_proc!(diags_proc, proc, ac, x, y, z, ::CPUs)
    for i in 1:size(ac,1)
        @inbounds diags_proc[x[i], y[i], z[i]] += proc[i] * ac[i]
    end
end

function diags_spcs!(diags_sp, proc, plank, ac, x, y, z, arch::Architecture)
    for diag in keys(diags_sp)
        if diag == :num || diag == :graz || diag == :mort || diag == :dvid
            nothing
        elseif diag == :PS
            diags_proc!(diags_sp[diag], proc.PS, ac, x, y, z, arch)
        elseif diag == :BS
            diags_proc!(diags_sp[diag], proc.BS, ac, x, y, z, arch)
        elseif diag == :VDOC
            diags_proc!(diags_sp[diag], proc.VDOC, ac, x, y, z, arch)
        elseif diag == :VNH4
            diags_proc!(diags_sp[diag], proc.VNH4, ac, x, y, z, arch)
        elseif diag == :VNO3
            diags_proc!(diags_sp[diag], proc.VNO3, ac, x, y, z, arch)
        elseif diag == :VPO4
            diags_proc!(diags_sp[diag], proc.VPO4, ac, x, y, z, arch)
        elseif diag == :resp
            diags_proc!(diags_sp[diag], proc.resp, ac, x, y, z, arch)
        elseif diag == :exu
            diags_proc!(diags_sp[diag], proc.exu, ac, x, y, z, arch)
        elseif diag == :Bm
            diags_proc!(diags_sp[diag], plank.Bm, ac, x, y, z, arch)
        elseif diag == :Cq
            diags_proc!(diags_sp[diag], plank.Cq, ac, x, y, z, arch)
        elseif diag == :Nq
            diags_proc!(diags_sp[diag], plank.Nq, ac, x, y, z, arch)
        elseif diag == :Pq
            diags_proc!(diags_sp[diag], plank.Pq, ac, x, y, z, arch)
        elseif diag == :chl
            diags_proc!(diags_sp[diag], plank.chl, ac, x, y, z, arch)
        else
            throw(ArgumentError("$(diag) is not one of the diagnostics"))
        end
    end
end
