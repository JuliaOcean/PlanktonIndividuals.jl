"""
    nutrients_init(g)
"""
function nutrients_init(::CPUs,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    nut = nutrient_fields(zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size),
                          zeros(total_size), zeros(total_size))
    return nut
end
function nutrients_init(::GPUs,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    nut = nutrient_fields(CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)),
                          CuArray(zeros(total_size)), CuArray(zeros(total_size)))
    return nut
end

"""
    gen_nutrients(g,nut)
Set up initial nutrient fields according to grid information
'nut' is an array of 10 elements, each element is a kind of nutrient
"""
function gen_nutrients(::CPUs,g,nut)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    DIC = fill(nut[1],total_size) .* rand(Uniform(0.8,1.2), total_size)
    NH4 = fill(nut[2],total_size) .* rand(Uniform(0.8,1.2), total_size)
    NO3 = fill(nut[3],total_size) .* rand(Uniform(0.8,1.2), total_size)
    PO4 = fill(nut[4],total_size) .* rand(Uniform(0.8,1.2), total_size)
    DOC = fill(nut[5],total_size) .* rand(Uniform(0.8,1.2), total_size)
    DON = fill(nut[6],total_size) .* rand(Uniform(0.8,1.2), total_size)
    DOP = fill(nut[7],total_size) .* rand(Uniform(0.8,1.2), total_size)
    POC = fill(nut[8],total_size) .* rand(Uniform(0.8,1.2), total_size)
    PON = fill(nut[9],total_size) .* rand(Uniform(0.8,1.2), total_size)
    POP = fill(nut[10],total_size) .* rand(Uniform(0.8,1.2), total_size)

    nutrients = nutrient_fields(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP)
    fill_halo!(nutrients,g)
    return nutrients
end
function gen_nutrients(::GPUs,g,nut)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    DIC = fill(nut[1],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    NH4 = fill(nut[2],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    NO3 = fill(nut[3],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    PO4 = fill(nut[4],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    DOC = fill(nut[5],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    DON = fill(nut[6],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    DOP = fill(nut[7],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    POC = fill(nut[8],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    PON = fill(nut[9],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray
    POP = fill(nut[10],total_size) .* rand(Uniform(0.8,1.2), total_size) |> CuArray

    nutrients = nutrient_fields(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP)
    fill_halo!(nutrients,g)
    return nutrients
end

"""
    load_nut_initials(paths,g)
Load nutrient initial conditions from files
"""
function load_nut_initials(::CPUs,paths,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    indices = ["DIC", "NH4", "NO3", "PO4", "DOC", "DON", "DOP", "POC", "PON", "POP"]
    pathkeys = collect(keys(paths))
    tmps = []
    for index in indices
        if length(findall(x->x==index, pathkeys)) == 0
            print("NUT_INIT: nutrient not found \n")
        else
            tmp = deserialize(paths[index])
            if size(tmp) == (g.Nx, g.Ny, g.Nz)
                tt = zeros(total_size)
                tt[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= tmp
                push!(tmps,tt)
            else
                print("NUT_INIT: grid mismatch \n")
            end
        end
    end
    nut = nutrient_fields(tmps[1],tmps[2],tmps[3],tmps[4],tmps[5],tmps[6],tmps[7],tmps[8],tmps[9],tmps[10])
    fill_halo!(nut,g)
    return nut
end
function load_nut_initials(::GPUs,paths,g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    indices = ["DIC", "NH4", "NO3", "PO4", "DOC", "DON", "DOP", "POC", "PON", "POP"]
    pathkeys = collect(keys(paths))
    tmps = []
    for index in indices
        if length(findall(x->x==index, pathkeys)) == 0
            print("NUT_INIT: nutrient not found \n")
        else
            tmp = deserialize(paths[index]) |> CuArray
            if size(tmp) == total_size
                tt = zeros(total_size) |> CuArray
                tt[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] .= tmp
                push!(tmps,tmp)
            else
                print("NUT_INIT: grid mismatch \n")
            end
        end
    end
    nut = nutrient_fields(tmps[1],tmps[2],tmps[3],tmps[4],tmps[5],tmps[6],tmps[7],tmps[8],tmps[9],tmps[10])
    fill_halo!(nut,g)
    return nut
end
