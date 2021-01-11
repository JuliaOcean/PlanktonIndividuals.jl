##### deactivate grazed or dead individuals
function deactivate!(plank, loss)
    @inbounds plank.ac .*= (1.0 .- loss)
end

##### grazing and grazing loss
function grazing!(plank, arch::Architecture, plk, p)
    ##### calculate grazing loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.graz, 
               p.grazFracC, p.grazFracN, p.grazFracP, p.R_NC, p.R_PC, arch)
    
    ##### deactivate grazed individuals
    deactivate!(plank, plank.graz)

    return nothing
end

##### mortality and mortality loss
function mortality!(plank, arch::Architecture, plk, p)
    ##### calculate mortality loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.mort, 
               p.mortFracC, p.mortFracN, p.mortFracP, p.R_NC, p.R_PC, arch)
    
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
    event = kernel!(idx, con, con_ind, de_ind)
    wait(device(arch), event)
    return nothing
end

##### copy ready to divide individuals to inactive rows
@kernel function copy_daughter_individuals_kernel!(plank, con, idx)
    i = @index(Global, Linear)
    if (con[i] == 1.0) & (idx[i] ≠ 0)
        # @print("index: $(idx[i]), $i \n")
        @inbounds plank.x[idx[i]]    = plank.x[i]
        @inbounds plank.y[idx[i]]    = plank.y[i]
        @inbounds plank.z[idx[i]]    = plank.z[i]
        @inbounds plank.Sz[idx[i]]   = plank.Sz[i]
        @inbounds plank.Bm[idx[i]]   = plank.Bm[i]
        @inbounds plank.Cq[idx[i]]   = plank.Cq[i]
        @inbounds plank.Nq[idx[i]]   = plank.Nq[i]
        @inbounds plank.Pq[idx[i]]   = plank.Pq[i]
        @inbounds plank.chl[idx[i]]  = plank.chl[i]
        @inbounds plank.gen[idx[i]]  = plank.gen[i]
        @inbounds plank.ac[idx[i]]   = plank.ac[i]
        @inbounds plank.dvid[idx[i]] = plank.dvid[i]
        @inbounds plank.graz[idx[i]] = plank.graz[i]
        @inbounds plank.mort[idx[i]] = plank.mort[i]
    end
end
function copy_daughter_individuals!(plank, con, idx::AbstractArray{Int64,1}, arch)
    kernel! = copy_daughter_individuals_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, con, idx)
    wait(device(arch), event)
    return nothing
end

##### cell division
@kernel function divide_to_half_kernel!(plank)
    i = @index(Global)
    @inbounds plank.Sz[i]  *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.Bm[i]  *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.Cq[i]  *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.Nq[i]  *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.Pq[i]  *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.chl[i] *= (2.0 - plank.dvid[i]) / 2 
    @inbounds plank.gen[i] += plank.dvid[i]
    @inbounds plank.age[i] *= (1.0 - plank.dvid[i])
    @inbounds plank.iS[i]   = plank.iS[i] * (1.0 - plank.dvid[i]) + plank.Sz[i] * plank.dvid[i]
end
function divide_to_half!(plank, arch)
    kernel! = divide_to_half_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank)
    wait(device(arch), event)
    return nothing
end
function divide!(plank, deactive_ind, arch::Architecture)
    con_ind = cumsum(plank.dvid)
    get_tind!(plank.idx, plank.dvid, Int.(con_ind), deactive_ind, arch)
    copy_daughter_individuals!(plank, plank.dvid, Int.(plank.idx), arch)
    divide_to_half!(plank, arch)
    # println("time step ends")
    return nothing
end