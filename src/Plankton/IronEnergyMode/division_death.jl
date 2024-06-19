##### inactivate grazed or dead individuals
function inactivate!(plank, loss)
    @inbounds plank.ac .*= (1.0f0 .- loss)
end

##### grazing and grazing loss
function grazing!(plank, arch::Architecture, plk, p)
    ##### calculate grazing loss
    calc_loss!(plk.NH4.data, plk.NO3.data, plk.DOC.data, plk.POC.data, 
               plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               plk.DOFe.data, plk.POFe.data, plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.graz, 
               p.grazFracC, p.grazFracN, p.grazFracP, p.grazFracFe, p.R_NC, p.R_PC, arch)
    
    ##### inactivate grazed individuals
    inactivate!(plank, plank.graz)

    return nothing
end

##### mortality and mortality loss
function mortality!(plank, arch::Architecture, plk, p)
    ##### calculate mortality loss
    calc_loss!(plk.NH4.data, plk.NO3.data, plk.DOC.data, plk.POC.data, 
               plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               plk.DOFe.data, plk.POFe.data, plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.graz, 
               p.mortFracC, p.mortFracN, p.mortFracP, p.mortFracFe, p.R_NC, p.R_PC, arch)
    
    ##### inactivate dead individuals
    inactivate!(plank, plank.mort)

    return nothing
end

@kernel function get_tind_kernel!(idx, con, con_ind, de_ind)
    i = @index(Global, Linear)
    if con[i] == 1.0f0
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
    if (con[i] == 1.0f0) & (idx[i] ≠ 0)
        # @print("index: $(idx[i]), $i \n")
        @inbounds plank.x[idx[i]]    = plank.x[i]
        @inbounds plank.y[idx[i]]    = plank.y[i]
        @inbounds plank.z[idx[i]]    = plank.z[i]
        @inbounds plank.Sz[idx[i]]   = plank.Sz[i]
        @inbounds plank.Bm[idx[i]]   = plank.Bm[i]
        @inbounds plank.CH[idx[i]]   = plank.CH[i]
        @inbounds plank.qNO3[idx[i]] = plank.qNO3[i]
        @inbounds plank.qNH4[idx[i]] = plank.qNH4[i]
        @inbounds plank.qP[idx[i]]   = plank.qP[i]
        @inbounds plank.qFe[idx[i]]  = plank.qFe[i]
        @inbounds plank.Chl[idx[i]]  = plank.Chl[i]
        @inbounds plank.gen[idx[i]]  = plank.gen[i]
        @inbounds plank.ac[idx[i]]   = plank.ac[i]
        @inbounds plank.dvid[idx[i]] = plank.dvid[i]
        @inbounds plank.graz[idx[i]] = plank.graz[i]
        @inbounds plank.mort[idx[i]] = plank.mort[i]
    end
end
function copy_daughter_individuals!(plank, con, idx::AbstractArray{Int,1}, arch)
    kernel! = copy_daughter_individuals_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, con, idx)
    return nothing
end

##### cell division
@kernel function divide_to_half_kernel!(plank)
    i = @index(Global)
    @inbounds plank.Sz[i]   *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.Bm[i]   *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.CH[i]   *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.qNO3[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.qNH4[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.qP[i]   *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.qFe[i]  *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.Chl[i]  *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.gen[i]  += plank.dvid[i]
    @inbounds plank.age[i]  *= (1.0f0 - plank.dvid[i])
end
function divide_to_half!(plank, arch)
    kernel! = divide_to_half_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank)
    return nothing
end
function divide!(plank, nuts, deactive_ind, arch::Architecture)
    accumulate!(+, nuts.idc, Int.(plank.dvid))
    get_tind!(plank.idx, plank.dvid, nuts.idc, deactive_ind, arch)
    copy_daughter_individuals!(plank, plank.dvid, plank.idx, arch)
    divide_to_half!(plank, arch)    
    return nothing
end
