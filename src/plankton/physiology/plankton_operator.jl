##### copy active individuals to model.timestepper.tmp
@kernel function copyto_tmp_kernel!(plank, tmp, con, idx, b::Bool)
    i = @index(Global, Linear)
    if con[i] == 1.0
        @inbounds tmp.x[idx[i]]    = copy(plank.x[i])
        @inbounds tmp.y[idx[i]]    = copy(plank.y[i])
        @inbounds tmp.z[idx[i]]    = copy(plank.z[i])
        @inbounds tmp.xi[idx[i]]   = copy(plank.xi[i])
        @inbounds tmp.yi[idx[i]]   = copy(plank.yi[i])
        @inbounds tmp.zi[idx[i]]   = copy(plank.zi[i])
        @inbounds tmp.iS[idx[i]]   = copy(plank.iS[i])
        @inbounds tmp.Sz[idx[i]]   = copy(plank.Sz[i])
        @inbounds tmp.Bm[idx[i]]   = copy(plank.Bm[i])
        @inbounds tmp.Cq[idx[i]]   = copy(plank.Cq[i])
        @inbounds tmp.Nq[idx[i]]   = copy(plank.Nq[i])
        @inbounds tmp.Pq[idx[i]]   = copy(plank.Pq[i])
        @inbounds tmp.chl[idx[i]]  = copy(plank.chl[i])
        @inbounds tmp.gen[idx[i]]  = copy(plank.gen[i])
        @inbounds tmp.age[idx[i]]  = copy(plank.age[i])
        @inbounds tmp.ac[idx[i]]   = copy(plank.ac[i])
        @inbounds tmp.graz[idx[i]] = copy(plank.graz[i])
        @inbounds tmp.mort[idx[i]] = copy(plank.mort[i])
        @inbounds tmp.dvid[idx[i]] = copy(plank.dvid[i])

        ##### (de)activate individuals
        @inbounds plank.ac[i] = plank.ac[i] * b
    end
end
function copyto_tmp!(plank, tmp, con, idx::AbstractArray{Int64,1}, b::Bool, arch)
    kernel! = copyto_tmp_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, tmp, con, idx, b)
    wait(device(arch), event)
    return nothing
end
##### grazing and grazing loss
function grazing!(plank, tmp, arch::Architecture, g::Grids, plk, p)
    ##### calculate index for timestepper.tmp
    plank.idx .= cumsum(plank.graz)
    ##### copy grazed individuals to timestepper.tmp
    copyto_tmp!(plank, tmp, plank.graz, Int.(plank.idx), false, arch)
    ##### calculate grazing loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               tmp, Int.(tmp.xi), Int.(tmp.yi), Int.(tmp.zi), 
               p.grazFracC, p.grazFracN, p.grazFracP, p.R_NC, p.R_PC, g, arch)
end

##### mortality and mortality loss
function mortality!(plank, tmp, arch::Architecture, g::Grids, plk, p)
    ##### calculate index for timestepper.tmp
    plank.idx .= cumsum(plank.mort)
    ##### copy dead individuals to timestepper.tmp
    copyto_tmp!(plank, tmp, plank.mort, Int.(plank.idx), false, arch)
    ##### calculate mortality loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               tmp, Int.(tmp.xi), Int.(tmp.yi), Int.(tmp.zi), 
               p.mortFracC, p.mortFracN, p.mortFracP, p.R_NC, p.R_PC, g, arch)
end

##### cell division
@kernel function divide_half_kernel!(plank, con)
    i = @index(Global, Linear)
    if con[i] == 1.0
        @inbounds plank.Sz[i]  *= 0.45
        @inbounds plank.Bm[i]  *= 0.45
        @inbounds plank.Cq[i]  *= 0.5
        @inbounds plank.Nq[i]  *= 0.5
        @inbounds plank.Pq[i]  *= 0.5
        @inbounds plank.chl[i] *= 0.5
        @inbounds plank.gen[i] += 1.0
        @inbounds plank.age[i]  = 1.0
        @inbounds plank.iS[i]   = copy(plank.Sz[i])
    end
end
function divide_to_half!(plank, con, arch::Architecture)
    kernel! = divide_half_kernel!(device(arch), 256, (size(plank.ac,1),))
    event = kernel!(plank, con)
    wait(device(arch), event)
    return nothing
end

function divide!(plank, arch::Architecture)
    ##### perform cell division
    divide_to_half!(plank, plank.dvid, arch)

    ##### calculate index for plank.data
    plank.idx .= cumsum(plank.dvid)
    plank.idx .+= dot(plank.ac, plank.ac)

    ##### copy newly divided individuals to the end of active individuals in plank.data
    copyto_tmp!(plank, plank, plank.dvid, Int.(plank.idx), true, arch)
    plank.dvid .= 0.0
end

##### zero a StructArray of plank.data
function zero_tmp!(tmp)
    tmp.x   .= 0.0
    tmp.y   .= 0.0
    tmp.z   .= 0.0
    tmp.xi  .= 0.0
    tmp.yi  .= 0.0
    tmp.zi  .= 0.0
    tmp.iS  .= 0.0
    tmp.Sz  .= 0.0
    tmp.Bm  .= 0.0
    tmp.Cq  .= 0.0
    tmp.Nq  .= 0.0
    tmp.Pq  .= 0.0
    tmp.chl .= 0.0
    tmp.gen .= 0.0
    tmp.age .= 0.0
    tmp.ac  .= 0.0
    tmp.idx .= 0.0
    tmp.graz .= 0.0
    tmp.mort .= 0.0
    tmp.dvid .= 0.0
    return nothing
end
