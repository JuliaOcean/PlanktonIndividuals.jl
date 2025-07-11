function tracers_init(arch, g, FT = Float32)
    fields = (Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT))

    tracers = NamedTuple{tracer_names}(fields)
    return tracers
end

"""
    default_tracer_init()
Generate defalut bgc tracer initial conditions.
"""
function default_tracer_init()
    init = (DIC=20.0, NH4=0.5, NO3=0.8, PO4=0.10, FeT=1.0e-6, DOC=10.0, DON=0.1, DOP=0.05, DOFe=0.0, POC=0.0, PON=0.0, POP=0.0, POFe=0.0, CHO = 0.1)
    rand_noise = (DIC=0.1, NH4=0.1, NO3=0.1, PO4=0.1, FeT=0.1, DOC=0.1, DON=0.1, DOP=0.1, DOFe=0.1, POC=0.1, PON=0.1, POP=0.1, POFe=0.1, CHO = 0.1)
    return (initial_condition = init, rand_noise = rand_noise)
end

"""
    generate_tracers(arch, grid, source, FT)
Set up initial bgc tracer fields according to `grid`.

Arguments
=================
- `arch`: `CPU()` or `GPU()`. The computer architecture used to time-step `model`.
- `grid`: The resolution and discrete geometry on which nutrient fields are solved.
- `source`: A `NamedTuple` containing 10 numbers each of which is the uniform initial 
            condition of one tracer, or a `Dict` containing the file paths pointing to
            the files of nutrient initial conditions.
- `FT`: Floating point data type. Default: `Float32`.
"""
function generate_tracers(arch::Architecture, g::AbstractGrid, 
                            source::Union{Dict,NamedTuple}, FT::DataType)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    tracers = tracers_init(arch, g, FT)
    pathkeys = collect(keys(source))

    if typeof(source) <: NamedTuple
        pathkeys = collect(keys(source.initial_condition))
    end

    for name in tracer_names
        if length(findall(x->x==name, pathkeys)) == 0
            throw(ArgumentError("TRAC_INIT: tracer not found $(name)"))
        else
            if typeof(source) <: Dict # file paths
                tmp = deserialize(source[name]) |> array_type(arch)
                if sum(tmp .< 0.0) > 0
                    throw(ArgumentError("TRAC_INIT: The initial condition should be none-negetive."))
                end

                if size(tmp) == (g.Nx, g.Ny, g.Nz)
                    @views @. tracers[name].data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] = FT.(tmp[:,:,:])
                else
                    throw(ArgumentError("TRAC_INIT:  grid mismatch"))
                end
            elseif typeof(source) <: NamedTuple # value
                if source.initial_condition[name] < 0.0
                    throw(ArgumentError("TRAC_INIT:  The initial condition should be none-negetive."))
                end
                lower = FT(1.0 - source.rand_noise[name])
                upper = FT(1.0 + source.rand_noise[name])
                tracers[name].data .= fill(FT(source.initial_condition[name]),total_size) .* rand(lower:1.0f-4:upper, total_size) |> array_type(arch)
            end
        end

        @views @. tracers[name].data *= g.landmask
    end

    fill_halo_tracer!(tracers,g)

    return tracers
end
