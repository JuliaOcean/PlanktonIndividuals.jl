##### compounds like exo-enzymes/signaling compounds will be adsorbed to abiotic particles
@kernel function calc_adsorption_kernel!(abiotics, cho, p)
  i = @index(Global)
  @inbounds abiotics.CHOe[i] = max(0.0f0, cho[abiotics.xi[i],abiotics.yi[i],abiotics.zi[i]]) * abiotics.ac[i]
  @inbounds abiotics.ADS[i] = abiotics.CHOe[i] * p.k_ads
end
function calc_adsorption!(abiotics, cho, p, arch)
  kernel! = calc_adsorption_kernel!(device(arch), 256, (size(abiotics.ac,1)))
  kernel!(abiotics, cho, p)
  return nothing
end

##### compounds might decay after adsorpted to abiotic particles
@kernel function calc_decay_kernel!(abiotics, p)
  i = @index(Global)
  @inbounds abiotics.DEC[i] = abiotics.CHO[i] * p.k_decay
end
function calc_decay!(abiotics, p, arch)
  kernel! = calc_decay_kernel!(device(arch), 256, (size(abiotics.ac,1)))
  kernel!(abiotics, p)
  return nothing
end

##### update C quotas
@kernel function update_quotas_kernel!(abiotics, ΔT)
  i = @index(Global)
  @inbounds abiotics.CHO[i]  += (abiotics.ADS[i] - abiotics.DEC[i]) * ΔT
end
function update_quotas!(abiotics, ΔT, arch)
  kernel! = update_quotas_kernel!(device(arch), 256, (size(abiotics.ac,1)))
  kernel!(abiotics, ΔT)
  return nothing
end

##### deal with chemical compounds adsorption
@kernel function calc_consume_kernel!(ctscho, abiotics, ac, x, y, z, ΔT)
  i = @index(Global)
  @inbounds KernelAbstractions.@atomic ctscho[x[i], y[i], z[i]] += -abiotics.ADS[i] * ΔT * ac[i]
end
function calc_consume!(ctscho, abiotics, ac, x, y, z, ΔT, arch)
  kernel! = calc_consume_kernel!(device(arch), 256, (size(ac,1)))
  kernel!(ctscho, abiotics, ac, x, y, z, ΔT)
  return nothing 
end

function abiotic_particle_update!(abiotic, cho, ctscho, diags_spcs, ΔT, arch::Architecture)
  particle = abiotic.data
  p = abiotic.p
  calc_adsorption!(particle, cho, p, arch)
  calc_decay!(particle, p, arch)
  update_quotas!(particle, ΔT, arch)
  calc_consume!(ctscho, particle, particle.ac, particle.xi, particle.yi, particle.zi, ΔT, arch)

  ##### Diagnostics 
  diags_spcs!(diags_spcs, abiotic, particle.ac, particle.xi, particle.yi, particle.zi, arch)
  ##### diagnostic for individual distribution
  diags_proc!(diags_spcs.num, particle.ac, particle.ac, particle.xi, particle.yi, particle.zi, arch)
end