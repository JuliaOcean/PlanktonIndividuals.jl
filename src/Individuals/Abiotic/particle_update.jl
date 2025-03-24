##### compounds like anti-biotics will be adsorbed to abiotic particles
@kernel function calc_adsorption_kernel!(particle, cho, p)
  i = @index(Global)
  @inbounds particle.CHOe[i] = max(0.0f0, cho[particle.xi[i],particle.yi[i],particle.zi[i]]) * particle.ac[i]
  @inbounds particle.ADS[i] = particle.CHOe[i] * p.k_ads
end
function calc_adsorption!(particle, cho, p, arch)
  kernel! = calc_adsorption_kernel!(device(arch), 256, (size(particle.ac,1)))
  kernel!(particle, cho, p)
  return nothing
end

##### compounds might decay after adsorpted to abiotic particles
@kernel function calc_decay_kernel!(particle, p)
  i = @index(Global)
  @inbounds particle.DEC[i] = particle.CHO[i] * p.k_decay
end
function calc_decay!(particle, p, arch)
  kernel! = calc_decay_kernel!(device(arch), 256, (size(particle.ac,1)))
  kernel!(particle, p)
  return nothing
end

##### update C quotas
@kernel function update_quotas_kernel!(particle, ΔT)
  i = @index(Global)
  @inbounds particle.CHO[i]  += (particle.ADS[i] - particle.DEC[i]) * ΔT
end
function update_quotas!(particle, ΔT, arch)
  kernel! = update_quotas_kernel!(device(arch), 256, (size(particle.ac,1)))
  kernel!(particle, ΔT)
  return nothing
end

##### deal with chemical compounds adsorption
@kernel function calc_consume_kernel!(ctscho, particle, ac, x, y, z, ΔT)
  i = @index(Global)
  @inbounds KernelAbstractions.@atomic ctscho[x[i], y[i], z[i]] += -particle.ADS[i] * ΔT * ac[i]
end
function calc_consume!(ctscho, particle, ac, x, y, z, ΔT, arch)
  kernel! = calc_consume_kernel!(device(arch), 256, (size(ac,1)))
  kernel!(ctscho, particle, ac, x, y, z, ΔT)
  return nothing 
end

function particle_update!(particle, cho, ctscho, p, diags_spcs, ΔT, arch::Architecture)
  calc_adsorption!(particle, cho, p, arch)
  calc_decay!(particle, p, arch)
  update_quotas!(particle, ΔT, arch)
  calc_consume!(ctscho, particle, particle.ac, particle.xi, particle.yi, particle.zi, ΔT, arch)

  ##### Diagnostics

end