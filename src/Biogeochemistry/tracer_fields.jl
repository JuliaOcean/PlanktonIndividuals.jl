struct Field{FT}
    data::AbstractArray{FT,3}
    bc::BoundaryConditions
end
"""
    Field(arch::Architecture, grid::AbstractGrid, FT::DataType; bcs = default_bcs())
Construct a `Field` on `grid` with data and boundary conditions on architecture `arch`
with DataType `FT`.
"""
function Field(arch::Architecture, grid::AbstractGrid, FT::DataType; bcs = default_bcs())
    total_size = (grid.Nx+grid.Hx*2, grid.Ny+grid.Hy*2, grid.Nz+grid.Hz*2)
    data = zeros(FT, total_size) |> array_type(arch)
    return Field{FT}(data,bcs)
end

@inline interior(c, grid) = c[grid.Hx+1:grid.Hx+grid.Nx, grid.Hy+1:grid.Hy+grid.Ny, grid.Hz+1:grid.Hz+grid.Nz]

function zero_fields!(a)
    for tr in keys(a)
        @inbounds a[tr].data .= 0.0f0
    end
end

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
    init = (DIC=20.0, NH4=0.5, NO3=0.8, PO4=0.10, DFe=1.0e-6, DOC=10.0, DON=0.1, DOP=0.05, PFe_inorg=0.0, POC=0.0, PON=0.0, POP=0.0, PFe_bio=0.0, Dust=0.0)
    rand_noise = (DIC=0.1, NH4=0.1, NO3=0.1, PO4=0.1, DFe=0.1, DOC=0.1, DON=0.1, DOP=0.1, PFe_inorg=0.1, POC=0.1, PON=0.1, POP=0.1, PFe_bio=0.1, Dust=0.1)
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
