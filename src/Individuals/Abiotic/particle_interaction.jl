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
@kernel function calc_merge_kernel!(merge_matrix, abiotic, phyto, grid, p)
  i, j = @index(Global, NTuple)
  merge_matrix[i,j] = isless(calc_distance(phyto.xi[i], phyto.yi[i],
                                           phyto.zi[i], phyto.x[i],
                                           phyto.y[i], phyto.z[i],
                                           abiotic.xi[j], abiotic.yi[j],
                                           abiotic.zi[j], abiotic.x[j],
                                           abiotic.y[j], abiotic.z[j],
                                           grid), p.Rd) * abiotic.ac[j] 
end
function calc_merge!(merge_matrix, abiotic, phyto, grid, p, arch::Architecture)
  kernel! = calc_merge_kernel!(device(arch), (16,16), (size(merge_matrix)))
  kernel!(merge_matrix, abiotic, phyto, grid, p)
  return nothing
end

##### clean the merge matrix - one vesicle might be close enough
##### to two or more phytoplankton cells.
##### Only the first phytoplankton cell can uptake this vesicle.


##### phytoplankton cell uptake the content of vesicles
@kernel function calc_CHOuptake_kernel!(abiotic, phytoCHO)
  i = @index(Global)
  phytoCHO += abiotic.CHO[i] * abiotic.merg[i]
end
function calc_CHOuptake!(abiotic, phytoCHO, arch::Architecture)
  kernel! = calc_CHOuptake_kernel!(device(arch), 256, (size(abiotic.ac,1)))
  kernel!(abiotic, phytoCHO)
  return nothing
end

##### inactivate merged vesicles
function inactivate!(abiotic)
  @inbounds abiotic.ac .*= (1.0f0 .- abiotic.merg)
end

function abiotic_particle_merge!(abiotic, phyto, arch::Architecture)
  
  
end