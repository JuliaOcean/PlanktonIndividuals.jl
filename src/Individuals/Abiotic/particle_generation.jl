function construct_abiotic_particle(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x   = zeros(FT, maxN), y   = zeros(FT, maxN), z   = zeros(FT, maxN),
                          xi  = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN),
                          ac  = zeros(FT, maxN), merg= zeros(Int,maxN), idx = zeros(Int,maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :Rd, :release_P, :sz_min, :Ktr)

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

function initialize_abiotic_particle!(particle, N::Int, g::AbstractGrid, arch::Architecture)
    particle.data.ac[1:N] .= 1.0f0                                                                             # activity

    rand!(rng_type(arch), particle.data.x)
    rand!(rng_type(arch), particle.data.y)
    rand!(rng_type(arch), particle.data.z)

    particle.data.x  .=(particle.data.x .* g.Nx) .* particle.data.ac                                         # x, unit: grid spacing, starting from 0
    particle.data.y  .=(particle.data.y .* g.Ny) .* particle.data.ac                                         # y, unit: grid spacing, starting from 0
    particle.data.z  .=(particle.data.z .* g.Nz) .* particle.data.ac                                         # z, unit: grid spacing, starting from 0

    mask_individuals!(particle.data, g, N, arch)
end
