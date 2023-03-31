##### deactivate grazed or dead individuals
function deactivate!(plank, loss)
    @inbounds plank.ac .*= (1.0 .- loss)
end

##### grazing and grazing loss
function grazing!(plank, arch::Architecture, plk, p)
    ##### calculate grazing loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.graz, 
               p.grazFracC, p.grazFracN, p.grazFracP, p, arch)
    
    ##### deactivate grazed individuals
    deactivate!(plank, plank.graz)

    return nothing
end

##### mortality and mortality loss
function mortality!(plank, arch::Architecture, plk, p)
    ##### calculate mortality loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.mort, 
               p.mortFracC, p.mortFracN, p.mortFracP, p, arch)
    
    ##### deactivate dead individuals
    deactivate!(plank, plank.mort)

    return nothing
end

@kernel function get_tind_kernel!(idx, con, con_ind, de_ind)
    i = @index(Global, Linear)
    if con[i] == 1.0
        idx[i] = de_ind[con_ind[i]]
    end
end
function get_tind!(idx, con, con_ind, de_ind, arch)
    kernel! = get_tind_kernel!(device(arch), 256, (size(idx,1)))
    kernel!(idx, con, con_ind, de_ind)
    return nothing
end

##### copy ready to divide individuals to inactive rows
@kernel function copy_daughter_individuals_kernel!(plank, con, idx)
    i = @index(Global, Linear)
    if (con[i] == 1.0) & (idx[i] ≠ 0)
        @inbounds plank.x[idx[i]]    = plank.x[i]
        @inbounds plank.y[idx[i]]    = plank.y[i]
        @inbounds plank.z[idx[i]]    = plank.z[i]
        @inbounds plank.PRO[idx[i]]  = plank.PRO[i]
        @inbounds plank.DNA[idx[i]]  = plank.DNA[i]
        @inbounds plank.RNA[idx[i]]  = plank.RNA[i]
        @inbounds plank.CH[idx[i]]   = plank.CH[i]
        @inbounds plank.NST[idx[i]]  = plank.NST[i]
        @inbounds plank.PST[idx[i]]  = plank.PST[i]
        @inbounds plank.Chl[idx[i]]  = plank.Chl[i]
        @inbounds plank.gen[idx[i]]  = plank.gen[i]
        @inbounds plank.ac[idx[i]]   = plank.ac[i]
        @inbounds plank.dvid[idx[i]] = plank.dvid[i]
        @inbounds plank.graz[idx[i]] = plank.graz[i]
        @inbounds plank.mort[idx[i]] = plank.mort[i]
    end
end
function copy_daughter_individuals!(plank, con, idx::AbstractArray{Int64,1}, arch)
    kernel! = copy_daughter_individuals_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, con, idx)
    return nothing
end

##### cell division
@kernel function divide_to_half_kernel!(plank)
    i = @index(Global)
    @inbounds plank.PRO[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.DNA[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.RNA[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.CH[i]  *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.NST[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.PST[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.Chl[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.gen[i] += plank.dvid[i]
    @inbounds plank.age[i] *= (1.0 - plank.dvid[i])
end
function divide_to_half!(plank, arch)
    kernel! = divide_to_half_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank)
    return nothing
end
function divide!(plank, deactive_ind, arch::Architecture)
    con_ind = cumsum(plank.dvid)
    get_tind!(plank.idx, plank.dvid, Int.(con_ind), deactive_ind, arch)
    copy_daughter_individuals!(plank, plank.dvid, Int.(plank.idx), arch)
    divide_to_half!(plank, arch)
    return nothing
end