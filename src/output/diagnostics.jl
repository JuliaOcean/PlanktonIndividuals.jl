using KernelAbstractions.Extras.LoopInfo: @unroll

mutable struct Diagnostics
    spcs::NamedTuple       # for each species
    tr::NamedTuple         # for tracers
end

function diags_setup(ntrs, nprocs, g::Grids, Nsp::Int64, arch::Architecture)
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
@kernel function diags_num_kernel!(diags_num, ac, x, y, z, g::Grids)
    @unroll for i in 1:size(ac,1)
        @inbounds diags_num[x[i]+g.Hx, y[i]+g.Hy, z[i]+g.Hz] += 1.0 * ac[i]
    end
end
function diags_num!(diags_num, ac, x, y, z, g::Grids, arch::Architecture)
    kernel! = diags_num_kernel!(device(arch), 1, (1,))
    event = kernel!(diags_num, ac, x, y, z, g)
    wait(device(arch), event)
    return nothing
end

@kernel function diags_graz_kernel!(diags_graz, graz, ac, x, y, z, g::Grids)
    @unroll for i in 1:size(ac,1)
        @inbounds diags_graz[x[i]+g.Hx, y[i]+g.Hy, z[i]+g.Hz] += graz[i] * ac[i]
    end
end
function diags_graz!(diags_graz, graz, ac, x, y, z, g::Grids, arch::Architecture)
    kernel! = diags_graz_kernel!(device(arch), 1, (1,))
    event = kernel!(diags_graz, graz, ac, x, y, z, g)
    wait(device(arch), event)
    return nothing
end

@kernel function diags_mort_kernel!(diags_mort, mort, ac, x, y, z, g::Grids)
    @unroll for i in 1:size(ac,1)
        @inbounds diags_mort[x[i]+g.Hx, y[i]+g.Hy, z[i]+g.Hz] += mort[i] * ac[i]
    end
end
function diags_mort!(diags_mort, mort, ac, x, y, z, g::Grids, arch::Architecture)
    kernel! = diags_mort_kernel!(device(arch), 1, (1,))
    event = kernel!(diags_mort, mort, ac, x, y, z, g)
    wait(device(arch), event)
    return nothing
end

@kernel function diags_dvid_kernel!(diags_dvid, dvid, ac, x, y, z, g::Grids)
    @unroll for i in 1:size(ac,1)
        @inbounds diags_dvid[x[i]+g.Hx, y[i]+g.Hy, z[i]+g.Hz] += dvid[i] * ac[i]
    end
end
function diags_dvid!(diags_dvid, dvid, ac, x, y, z, g::Grids, arch::Architecture)
    kernel! = diags_dvid_kernel!(device(arch), 1, (1,))
    event = kernel!(diags_dvid, dvid, ac, x, y, z, g)
    wait(device(arch), event)
    return nothing
end

@kernel function diags_proc_kernel!(diags_proc, proc, ac, x, y, z, g::Grids)
    @unroll for i in 1:size(ac,1)
        @inbounds diags_proc[x[i]+g.Hx, y[i]+g.Hy, z[i]+g.Hz] += proc[i] * ac[i]
    end
end
function diags_proc!(diags_proc, proc, ac, x, y, z, g::Grids, arch::Architecture)
    kernel! = diags_proc_kernel!(device(arch), 1, (1,))
    event = kernel!(diags_proc, proc, ac, x, y, z, g)
    wait(device(arch), event)
    return nothing
end

function diags_spcs!(diags_sp, proc, plank, ac, x, y, z, g::Grids, arch::Architecture)
    for diag in keys(diags_sp)
        if diag == :PS
            diags_proc!(diags_sp[diag], proc.PS, ac, x, y, z, g, arch)
        elseif diag == :BS
            diags_proc!(diags_sp[diag], proc.BS, ac, x, y, z, g, arch)
        elseif diag == :VDOC
            diags_proc!(diags_sp[diag], proc.VODC, ac, x, y, z, g, arch)
        elseif diag == :VNH4
            diags_proc!(diags_sp[diag], proc.VNH4, ac, x, y, z, g, arch)
        elseif diag == :VNO3
            diags_proc!(diags_sp[diag], proc.VNO3, ac, x, y, z, g, arch)
        elseif diag == :VPO4
            diags_proc!(diags_sp[diag], proc.VPO4, ac, x, y, z, g, arch)
        elseif diag == :resp
            diags_proc!(diags_sp[diag], proc.resp, ac, x, y, z, g, arch)
        elseif diag == :exu
            diags_proc!(diags_sp[diag], proc.exu, ac, x, y, z, g, arch)
        elseif diag == :Bm
            diags_proc!(diags_sp[diag], plank.Bm, ac, x, y, z, g, arch)
        elseif diag == :Cq
            diags_proc!(diags_sp[diag], plank.Cq, ac, x, y, z, g, arch)
        elseif diag == :Nq
            diags_proc!(diags_sp[diag], plank.Nq, ac, x, y, z, g, arch)
        elseif diag == :Pq
            diags_proc!(diags_sp[diag], plank.Pq, ac, x, y, z, g, arch)
        elseif diag == :chl
            diags_proc!(diags_sp[diag], plank.chl, ac, x, y, z, g, arch)
        end
    end
end
