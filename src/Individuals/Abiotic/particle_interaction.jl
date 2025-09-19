##### calculate distance between particles only when the two
##### particles are in the same grid
@inline function calc_distance(xi1, yi1, zi1, x1, y1, z1,
                               xi2, yi2, zi2, x2, y2, z2, grid)
    if (xi1 == xi2) && (yi1 == yi2) && (zi1 == zi2)
        dx = (x1 - x2) * ΔxF(xi1, yi1, zi1, grid)
        dy = (y1 - y2) * ΔyF(xi1, yi1, zi1, grid)
        dz = (z1 - z2) * ΔzF(xi1, yi1, zi1, grid)
        dist = sqrt(dx^2.0f0 + dy^2.0f0 + dz^2.0f0)
   else
        dist = 1.0f20
    end
    return dist
end

##### compare particle distance and Rd, decide whether
##### this vesicle will be merged into phytoplankton cell
@kernel function calc_merge_matrix_kernel!(intac, abiotic, plank, grid, abio_p)
    i, j = @index(Global, NTuple)
    intac[i,j] = isless(calc_distance(plank.xi[i], plank.yi[i],
                                      plank.zi[i], plank.x[i],
                                      plank.y[i], plank.z[i],
                                      abiotic.xi[j], abiotic.yi[j],
                                      abiotic.zi[j], abiotic.x[j],
                                      abiotic.y[j], abiotic.z[j],
                                      grid), abio_p.Rd) * abiotic.ac[j] 
end
function calc_merge_matrix!(intac, abiotic, plank, grid, abio_p, arch::Architecture)
    kernel! = calc_merge_matrix_kernel!(device(arch), (16,16), (size(intac)))
    kernel!(intac, abiotic, plank, grid, abio_p)
    return nothing
end

##### clean the merge matrix - one vesicle might be close enough
##### to two or more phytoplankton cells. 
##### Only the first phytoplankton cell can uptake this vesicle.
##### First, find all indices of 1.0f0 in intac, return an AbstracArray of
##### Cartesian Index, e.g., CartesianIndex(phyto_ind, abiotic_ind)
##### Then, assign phyto_ind to abiotic.merge[abiotic_ind]
@kernel function assign_merge_kernel!(abiotic, inds)
    i = @index(Global)
    abiotic.merg[inds[i][2]] = inds[i][1]
end
function assign_merge!(abiotic, inds, arch::Architecture)
    kernel! = assign_merge_kernel!(device(arch), 256, (size(inds,1)))
    kernel!(abiotic, inds)
    return nothing
end

##### merge particles
@kernel function merge_particle_kernel!(abiotic, plank)
    i = @index(Global)
    if abiotic.merg[i] ≠ 0
        @inbounds KernelAbstractions.@atomic plank.ptc[abiotic.merg[i]] += 1.0f0
    end
end
function merge_particle!(abiotic, plank, arch::Architecture)
    kernel! = merge_particle_kernel!(device(arch), 256, (size(abiotic.ac,1)))
    kernel!(abiotic, plank)
    return nothing
end

##### particle interaction wrapper
function particle_interaction!(abiotic, plank, intac, abio_p, grid, arch::Architecture)
    calc_merge_matrix!(intac, abiotic, plank, grid, abio_p, arch)
    inds = findall(isequal(1.0f0), intac)
    assign_merge!(abiotic, inds, arch)
    merge_particle!(abiotic, plank, arch)
    unsafe_free!(inds)
    abiotic.ac .= isless.(abiotic.merg, 1)
    abiotic.merg .= 0
    return nothing
end

##### particle release from phytoplankton cells
function get_release_probability!(plank, rnd, abio_p, ΔT, arch::Architecture)
    ##### generate random number (0,1)
    rand!(rng_type(arch), rnd.x)

    ##### compare the random number with the given probability
    ##### return 1 if random number is smaller
    @inbounds plank.Rptc .= isless.(rnd.x, abio_p.release_P .* ΔT .* plank.ptc) .* plank.ac
    return nothing
end

@kernel function release_abiotic_particle_kernel!(plank, abiotic, con, idx)
    i = @index(Global)
    if (con[i] == 1.0f0) && (idx[i] ≠ 0)
        @inbounds abiotic.ac[idx[i]] = 1.0f0
        @inbounds abiotic.x[idx[i]]  = plank.x[i]
        @inbounds abiotic.y[idx[i]]  = plank.y[i]
        @inbounds abiotic.z[idx[i]]  = plank.z[i]
        @inbounds plank.ptc[i]      -= 1.0f0
    end
end
function release_abiotic_particle!(plank, abiotic, con, idx, arch::Architecture)
    kernel! = release_abiotic_particle_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, abiotic, con, idx)
    return nothing
end

function particle_release!(plank, abiotic, trs, rnd, abio_p, ΔT, t, arch::Architecture)
    get_release_probability!(plank, rnd, abio_p, ΔT, arch)
    releasenum = dot(plank.Rptc, plank.ac)
    deactive_ind = findall(x -> x == 0.0f0, abiotic.ac)
    if releasenum > length(deactive_ind)
        throw(ArgumentError("number of abiotic particles exceeds the capacity at timestep $(t/86400.0) days"))
    end
    accumulate!(+, trs.idc, plank.Rptc)
    trs.idc_int .= unsafe_trunc.(Int, trs.idc)
    get_tind!(plank.idx, plank.Rptc, trs.idc_int, deactive_ind, arch)
    release_abiotic_particle!(plank, abiotic, plank.Rptc, plank.idx, arch)
    plank.idx .= 0
    unsafe_free!(deactive_ind)
    return nothing
end

##### particle boundary conditions
##### only positive fluxes are accepted for particles
@kernel function calc_particle_bc_top_kernel!(tr_temp, abio_bc_top, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid)
    i, j = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = 1 + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min * abio_p.Nsuper, 
                                               getbc(abio_bc_top, i, j, iter) * ΔT * Az(ii, jj, kk, g))
end
@kernel function calc_particle_bc_bottom_kernel!(tr_temp, abio_bc_bottom, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid)
    i, j = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = g.Nz + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min * abio_p.Nsuper, 
                                               getbc(abio_bc_bottom, i, j, iter) * ΔT * Az(ii, jj, kk, g))
end
@kernel function calc_particle_bc_east_kernel!(tr_temp, abio_bc_east, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid)
    j, k = @index(Global, NTuple)
    ii = g.Nx + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min * abio_p.Nsuper, 
                                               getbc(abio_bc_east, j, k, iter) * ΔT * Ax(ii, jj, kk, g))
end
@kernel function calc_particle_bc_west_kernel!(tr_temp, abio_bc_west, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid)
    j, k = @index(Global, NTuple)
    ii = 1 + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min * abio_p.Nsuper, 
                                               getbc(abio_bc_west, j, k, iter) * ΔT * Ax(ii, jj, kk, g))
end
@kernel function calc_particle_bc_south_kernel!(tr_temp, abio_bc_south, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid)
    j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = 1 + g.Hy
    kk = k + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min * abio_p.Nsuper, 
                                               getbc(abio_bc_south, i, k, iter) * ΔT * Ay(ii, jj, kk, g))
end
@kernel function calc_particle_bc_north_kernel!(tr_temp, abio_bc_north, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid)
    j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = g.Ny + g.Hy
    kk = k + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min * abio_p.Nsuper, 
                                               getbc(abio_bc_north, i, k, iter) * ΔT * Ay(ii, jj, kk, g))
end

function calc_particle_bc_top!(tr_temp, abio_bc_top, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
    kernel! = calc_particle_bc_top_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    kernel!(tr_temp, abio_bc_top, rnd_3d, abio_p, ΔT, iter, g)
    return nothing
end
function calc_particle_bc_bottom!(tr_temp, abio_bc_bottom, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
    kernel! = calc_particle_bc_bottom_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    kernel!(tr_temp, abio_bc_bottom, rnd_3d, abio_p, ΔT, iter, g)
    return nothing
end
function calc_particle_bc_east!(tr_temp, abio_bc_east, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
    kernel! = calc_particle_bc_east_kernel!(device(arch), (16,16), (g.Ny, g.Nz))
    kernel!(tr_temp, abio_bc_east, rnd_3d, abio_p, ΔT, iter, g)
    return nothing
end
function calc_particle_bc_west!(tr_temp, abio_bc_west, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
    kernel! = calc_particle_bc_west_kernel!(device(arch), (16,16), (g.Ny, g.Nz))
    kernel!(tr_temp, abio_bc_west, rnd_3d, abio_p, ΔT, iter, g)
    return nothing
end
function calc_particle_bc_south!(tr_temp, abio_bc_south, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
    kernel! = calc_particle_bc_south_kernel!(device(arch), (16,16), (g.Nx, g.Nz))
    kernel!(tr_temp, abio_bc_south, rnd_3d, abio_p, ΔT, iter, g)
    return nothing
end
function calc_particle_bc_north!(tr_temp, abio_bc_north, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
    kernel! = calc_particle_bc_north_kernel!(device(arch), (16,16), (g.Nx, g.Nz))
    kernel!(tr_temp, abio_bc_north, rnd_3d, abio_p, ΔT, iter, g)
    return nothing
end

function calc_particle_bcs!(tr_temp, abio_bcs, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture)
       calc_particle_bc_top!(tr_temp, abio_bcs.top,    rnd_3d, abio_p, ΔT, iter, g, arch)
    calc_particle_bc_bottom!(tr_temp, abio_bcs.bottom, rnd_3d, abio_p, ΔT, iter, g, arch)
      calc_particle_bc_east!(tr_temp, abio_bcs.east,   rnd_3d, abio_p, ΔT, iter, g, arch)
      calc_particle_bc_west!(tr_temp, abio_bcs.west,   rnd_3d, abio_p, ΔT, iter, g, arch)
     calc_particle_bc_north!(tr_temp, abio_bcs.north,  rnd_3d, abio_p, ΔT, iter, g, arch)
     calc_particle_bc_south!(tr_temp, abio_bcs.south,  rnd_3d, abio_p, ΔT, iter, g, arch)
    return nothing
end
# Avoid some computation / memory accesses when no flux boundary conditions will be applied
   calc_particle_bc_top!(tr_temp, ::Nothing, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture) = nothing
calc_particle_bc_bottom!(tr_temp, ::Nothing, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture) = nothing
  calc_particle_bc_east!(tr_temp, ::Nothing, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture) = nothing
  calc_particle_bc_west!(tr_temp, ::Nothing, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture) = nothing
 calc_particle_bc_south!(tr_temp, ::Nothing, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture) = nothing
 calc_particle_bc_north!(tr_temp, ::Nothing, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, arch::Architecture) = nothing

@kernel function copy_abiotic_particle_from_field_kernel!(abiotic, inds, de_inds)
    i = @index(Global)
    @inbounds abiotic.ac[de_inds[i]] = 1.0f0
    @inbounds abiotic.x[de_inds[i]]  = inds[i][1] # add a random number (0,1) if necessary
    @inbounds abiotic.y[de_inds[i]]  = inds[i][2]
    @inbounds abiotic.z[de_inds[i]]  = inds[i][3]
end
function copy_abiotic_particle_from_field!(abiotic, inds, de_inds, arch::Architecture)
    kernel! = copy_abiotic_particle_from_field_kernel!(device(arch), 256, (size(inds,1)))
    kernel!(abiotic, inds, de_inds)
    return nothing
end

function particles_from_bcs!(abiotic, tr_temp, abio_bcs::BoundaryConditions, rnd_3d, abio_p, ΔT, iter, g::AbstractGrid, t, arch::Architecture)
    rand!(rng_type(arch), rnd_3d)
    tr_temp .= 0.0f0
    calc_particle_bcs!(tr_temp, abio_bcs, rnd_3d, abio_p, ΔT, iter, g, arch)
    if sum(tr_temp) == 0.0f0
        return nothing
    else
        inds = findall(isequal(1.0f0), tr_temp)
        deactive_ind = findall(x -> x == 0.0f0, abiotic.ac)
        if length(inds) > length(deactive_ind)
            throw(ArgumentError("number of abiotic particles exceeds the capacity at timestep $(t/86400.0) days"))
        end
        copy_abiotic_particle_from_field!(abiotic, inds, deactive_ind, arch)
        unsafe_free!(inds)
        unsafe_free!(deactive_ind)
        return nothing
    end
end

##### particle-field interaction
##### turn field tracer into particles
##### need to substract particles from the tracer field
#=
@kernel function tracer_to_particle_kernel!(tr_temp, tr, rnd_3d, abio_p, ΔT, g::AbstractGrid)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds tr_temp[ii, jj, kk] = isless(rnd_3d[ii, jj, kk] * abio_p.sz_min, 
                                           tr[ii, jj, kk] * abio_p.Ktr * ΔT * volume(ii, jj, kk, g))
end
function tracer_to_particle!(tr_temp, tr, rnd_3d, abio_p, ΔT, g::AbstractGrid, arch::Architecture)
    kernel! = tracer_to_particle_kernel!((device(arch), (16,16), (g.Nx, g.Ny, g.Nz)))
    kernel!(tr_temp, tr, rnd_3d, abio_p, ΔT, g)
    return nothing
end

function particles_from_tracer!(abiotic, tr_temp, tr, rnd_3d, abio_p, ΔT, g::AbstractGrid, arch::Architecture)
    if abio_p.Ktr == 0.0f0
        return nothing
    else
        rand!(rng_type(arch), rnd_3d)
        tracer_to_particle!(tr_temp, tr, rnd_3d, abio_p, ΔT, g, arch)
        inds = findall(isequal(1.0f0), tr_temp)
        deactive_ind = findall(x -> x == 0.0f0, abiotic.ac)
        if length(inds) > length(deactive_ind)
            throw(ArgumentError("number of abiotic particles exceeds the capacity at timestep $(t/86400.0) days"))
        end
        copy_abiotic_particle_from_tracer!(abiotic, inds, deactive_ind, arch)
        unsafe_free!(inds)
        unsafe_free!(deactive_ind)
        return nothing
    end
end

=#