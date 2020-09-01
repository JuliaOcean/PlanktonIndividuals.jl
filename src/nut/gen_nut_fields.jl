"""
    nutrients_init(arch, g)
generate `Fields` with `0.0` in `Field.data`
"""
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
    gen_nutrients(arch, g, nut)
Set up initial nutrient fields according to grid information
'nut' is an array of 10 elements, each element is a kind of nutrient
"""
function gen_nutrients(arch, g, nut)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    nutrients = nutrients_init(arch, g)

    for i in 1:size(nut_names)
        nutrients[i].data .= fill(nut[i],total_size) .* rand(Uniform(0.8,1.2), total_size) |> array_type(arch)
    end

    fill_halo!(nutrients,g)

    return nutrients
end

"""
    load_nut_initials(arch, paths, g)
Load nutrient initial conditions from files
"""
function load_nut_initials(arch, paths, g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)

    nut = nutrients_init(arch, g)

    pathkeys = collect(keys(paths))

    tmps = []

    for name in nut_names
        if length(findall(x->x==name, pathkeys)) == 0
            print("NUT_INIT: nutrient not found \n")
        else
            tmp = deserialize(paths[name])
            if size(tmp) == (g.Nx, g.Ny, g.Nz)
                nut[name].data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= tmp
            else
                print("NUT_INIT: grid mismatch \n")
            end
        end
    end

    fill_halo!(nut,g)
    return nut
end