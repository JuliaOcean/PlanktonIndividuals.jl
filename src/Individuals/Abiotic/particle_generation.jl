function construct_abiotic_particle(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x   = zeros(FT, maxN), y   = zeros(FT, maxN), z   = zeros(FT, maxN),
                          xi  = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN),
                          ac  = zeros(Bool, maxN), merg= zeros(Int,maxN), idx = zeros(Int,maxN),
                          Fe_con = zeros(FT, maxN), sz = zeros(FT, maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :Rd, :release_P, :sz_min, :Ktr, :ptc_de, :Fe_frac, :M_Fe)

    pkeys = Symbol.(collect(keys(params)))
    tmp = zeros(length(param_names))
    for i in 1:length(param_names)
        if param_names[i] ∉ pkeys
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(FT.(tmp))
    return abiotic_particle(data, p, default_bcs())
end


#### Calculate the radius based on the inorganic iron content.

@inline function calculate_iron_particles_radius(Fe_con::Float32, particle)
    mass_Fe = Fe_con * particle.M_Fe * 1e-6f0
    total_mass = mass_Fe / particle.Fe_frac
    volume = total_mass / particle.ptc_de
    radius = cbrt(volume * 3.0f0 / (4.0f0 * Float32(pi)))
    return radius
end

function initialize_abiotic_particle!(particle, N::Int, g::AbstractGrid, arch::Architecture)
    particle.data.ac[1:N] .= true       # activity

    particle.data.Fe_con[1:N] .= particle.p.sz_min
    particle.data.sz[1:N] .= calculate_iron_particles_radius.(particle.data.Fe_con[1:N], Ref(particle.p))

    rand!(rng_type(arch), particle.data.x)
    rand!(rng_type(arch), particle.data.y)
    rand!(rng_type(arch), particle.data.z)

    particle.data.x  .=(particle.data.x .* g.Nx) .* particle.data.ac                                         # x, unit: grid spacing, starting from 0
    particle.data.y  .=(particle.data.y .* g.Ny) .* particle.data.ac                                         # y, unit: grid spacing, starting from 0
    particle.data.z  .=(particle.data.z .* g.Nz) .* particle.data.ac                                         # z, unit: grid spacing, starting from 0

    mask_individuals!(particle.data, g, N, arch)
end
