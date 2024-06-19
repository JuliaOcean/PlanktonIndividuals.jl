function nutrients_init(arch, g, FT = Float32)
    fields = (Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT), Field(arch, g, FT),
              Field(arch, g, FT))

    nut = NamedTuple{nut_names}(fields)
    return nut
end

"""
    default_nut_init()
Generate defalut nutrient initial conditions.
"""
function default_nut_init()
    init = (DIC=20.0, NH4=0.5, NO3=0.8, PO4=0.10, FeT=1.0e-6, DOC=10.0, DON=0.1, DOP=0.05, DOFe=0.0, POC=0.0, PON=0.0, POP=0.0, POFe=0.0)
    rand_noise = (DIC=0.1, NH4=0.1, NO3=0.1, PO4=0.1, FeT=0.1, DOC=0.1, DON=0.1, DOP=0.1, DOFe=0.1, POC=0.1, PON=0.1, POP=0.1, POFe=0.1)
    return (initial_condition = init, rand_noise = rand_noise)
end

"""
    generate_nutrients(arch, grid, source, FT)
Set up initial nutrient fields according to `grid`.

Arguments
=================
- `arch`: `CPU()` or `GPU()`. The computer architecture used to time-step `model`.
- `grid`: The resolution and discrete geometry on which nutrient fields are solved.
- `source`: A `NamedTuple` containing 10 numbers each of which is the uniform initial 
            condition of one tracer, or a `Dict` containing the file paths pointing to
            the files of nutrient initial conditions.
- `FT`: Floating point data type. Default: `Float32`.
"""
function generate_nutrients(arch::Architecture, g::AbstractGrid, 
                            source::Union{Dict,NamedTuple}, FT::DataType)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    nut = nutrients_init(arch, g, FT)
    pathkeys = collect(keys(source))

    if typeof(source) <: NamedTuple
        pathkeys = collect(keys(source.initial_condition))
    end

    for name in nut_names
        if length(findall(x->x==name, pathkeys)) == 0
            throw(ArgumentError("NUT_INIT: nutrient not found $(name)"))
        else
            if typeof(source) <: Dict # file paths
                tmp = deserialize(source[name]) |> array_type(arch)
                if sum(tmp .< 0.0) > 0
                    throw(ArgumentError("NUT_INIT:  The initial condition should be none-negetive."))
                end

                if size(tmp) == (g.Nx, g.Ny, g.Nz)
                    @views @. nut[name].data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] = FT.(tmp[:,:,:])
                else
                    throw(ArgumentError("NUT_INIT:  grid mismatch"))
                end
            elseif typeof(source) <: NamedTuple # value
                if source.initial_condition[name] < 0.0
                    throw(ArgumentError("NUT_INIT:  The initial condition should be none-negetive."))
                end
                lower = FT(1.0 - source.rand_noise[name])
                upper = FT(1.0 + source.rand_noise[name])
                nut[name].data .= fill(FT(source.initial_condition[name]),total_size) .* rand(lower:1.0f-4:upper, total_size) |> array_type(arch)
            end
        end

        @views @. nut[name].data *= g.landmask
    end

    fill_halo_nut!(nut,g)

    return nut
end
