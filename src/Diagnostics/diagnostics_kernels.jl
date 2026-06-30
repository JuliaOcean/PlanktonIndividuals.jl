##### record diagnostics of particle processes at each time step
@kernel function diags_proc_kernel!(diags_proc, proc, ac, x, y, z)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic diags_proc[x[i], y[i], z[i]] += proc[i] * ac[i]
end
function diags_proc!(diags_proc, proc, ac, x, y, z, arch)
    kernel! = diags_proc_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(diags_proc, proc, ac, x, y, z)
    return nothing 
end

function diags_spcs!(diags_sp, plank::phytoplankton, ac, x, y, z, arch::Architecture)
    for diag in keys(diags_sp)
        if diag in (:num, :graz, :mort, :dvid, :ptc)
            nothing
        else
            diags_proc!(diags_sp[diag], getproperty(plank.data, diag), ac, x, y, z, arch)
        end
    end
end

function diags_colony!(diags_colony, colony, arch::Architecture)
    for sp in eachindex(colony)
        diags_spcs!(diags_colony[sp], colony[sp], colony[sp].data.ac, 
                    colony[sp].data.xi, colony[sp].data.yi, colony[sp].data.zi, arch)
    end
end