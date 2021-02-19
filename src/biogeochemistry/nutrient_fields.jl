function nutrients_init(arch, g)
    fields = (Field(arch, g), Field(arch, g),
              Field(arch, g), Field(arch, g),
              Field(arch, g), Field(arch, g),
              Field(arch, g), Field(arch, g),
              Field(arch, g), Field(arch, g))

    nut = NamedTuple{nut_names}(fields)
    return nut
end

"""
    generate_nutrients(arch, grid, source)
Set up initial nutrient fields according to `grid`.

Keyword Arguments
=================
- `arch`: `CPUs()` or `GPUs()`. The computer architecture used to time-step `model`.
- `grid`: The resolution and discrete geometry on which nutrient fields are solved.
- `source`: An 10-element array with each element representing the initial condition of a kind of nutrient, 
            or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
"""
function generate_nutrients(arch, g, source::Array)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    nutrients = nutrients_init(arch, g)

    for i in 1:length(nut_names)
        nutrients[i].data .= fill(source[i],total_size) .* rand(Uniform(0.8,1.2), total_size) |> array_type(arch)
    end

    fill_halo_nut!(nutrients,g)

    return nutrients
end

function generate_nutrients(arch, g, source::Dict)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)

    nut = nutrients_init(arch, g)

    pathkeys = collect(keys(source))

    tmps = []

    for name in nut_names
        if length(findall(x->x==name, pathkeys)) == 0
            throw(ArgumentError("NUT_INIT: nutrient not found $(name)"))
        else
            tmp = deserialize(source[name]) |> array_type(arch)
            if size(tmp) == (g.Nx, g.Ny, g.Nz)
                @views @. nut[name].data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] = tmp[:,:,:]
            else
                throw(ArgumentError("NUT_INIT:  grid mismatch"))
            end
        end
    end

    fill_halo_nut!(nut,g)
    return nut
end
