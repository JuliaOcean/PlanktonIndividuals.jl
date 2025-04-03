mutable struct PlanktonOutputWriter
    filepath::String
    write_log::Bool
    save_diags::Bool
    save_phytoplankton::Bool
    save_abiotic_particle::Bool
    diags_file::String
    phytoplankton_file::String
    phytoplankton_include::Tuple
    phytoplankton_iteration_interval::Int
    abiotic_particle_file::String
    abiotic_particle_include::Tuple
    abiotic_particle_iteration_interval::Int
    max_filesize::Number # in Bytes
    part_diags::Int
    part_phytoplankton::Int
    part_abiotic_particle::Int
end

"""
    PlanktonOutputWriter(;dir = "./results",
                               diags_prefix = "diags",
                               phytoplankton_prefix = "phytoplankton",
                               abiotic_particle_prefix = "abiotic_particle",
                               write_log = false,
                               save_diags = false,
                               save_phytoplankton = false,
                               save_abiotic_particle = false,
                               phytoplankton_include = (:x, :y, :z),
                               abiotic_particle_include = (:x, :y, :z),
                               phytoplankton_iteration_interval = 1,
                               abiotic_particle_iteration_interval = 1,
                               max_filesize = Inf,
                               )
Generate a `PlanktonOutputWriter` structure which includes settings for model outputs

Keyword Arguments (Optional)
============================
- `dir`: The directory to store model outputs, "./results" by default
- `diags_prefix`: Descriptive filename prefixed to diagnostic output files.
- `phytoplankton_prefix`: Descriptive filename prefixed to phytoplankton output files.
- `write_log`: write model logs which contain global averages of simulated phytoplankton, default: `false`.
- `save_diags`: write diagnostics to disk, default: `false`.
- `save_phytoplankton`: write phytoplankton to disk, default: `false`.
- `phytoplankton_include`: list of phytoplankton properties to save, default: `(:x, :y, :z, :Sz)`.
- `phytoplankton_iteration_interval`: The time interval that phytoplankton are saved, 1 timestep by default.
- `save_abiotic_particle`: write abiotic_particle to disk, default: `false`.
- `abiotic_particle_include`: list of abiotic_particle properties to save, default: `(:x, :y, :z, :Sz)`.
- `abiotic_particle_iteration_interval`: The time interval that abiotic_particle are saved, 1 timestep by default.
- `max_filesize`: The writer will stop writing to the output file once the file size exceeds `max_filesize`,
                    and write to a new one with a consistent naming scheme ending in `part1`, `part2`, etc.
                    default: `Inf`.
"""
function PlanktonOutputWriter(;dir = "./results",
                               diags_prefix = "diags",
                               phytoplankton_prefix = "phytoplankton",
                               abiotic_particle_prefix = "abiotic_particle",
                               write_log = false,
                               save_diags = false,
                               save_phytoplankton = false,
                               save_abiotic_particle = false,
                               phytoplankton_include = (:x, :y, :z),
                               abiotic_particle_include = (:x, :y, :z),
                               phytoplankton_iteration_interval = 1,
                               abiotic_particle_iteration_interval = 1,
                               max_filesize = Inf
                               )

    isdir(dir) && rm(dir, recursive=true)
    mkdir(dir)

    diags_file = ""
    phytoplankton_file = ""
    abiotic_particle_file = ""

    if save_diags
        diags_file = joinpath(dir, diags_prefix*".jld2")
    end
    if save_phytoplankton
        phytoplankton_file = joinpath(dir, phytoplankton_prefix*".jld2")
    end
    if save_abiotic_particle
        abiotic_particle_file = joinpath(dir, abiotic_particle_prefix*".jld2")
    end

    return PlanktonOutputWriter(dir, write_log, save_diags, save_phytoplankton,
                                save_abiotic_particle,
                                diags_file, phytoplankton_file, phytoplankton_include, 
                                phytoplankton_iteration_interval, abiotic_particle_file, 
                                abiotic_particle_include, abiotic_particle_iteration_interval, 
                                max_filesize, 1, 1, 1)
end

function show(io::IO, writer::PlanktonOutputWriter)
    print(io, "PlanktonOutputWriter:\n",
              "├── files are saved at $(writer.filepath)\n",
              "├── $(save_diags_string(writer))\n",
              "├── $(save_phyto_string(writer))\n",
              "├── $(save_abiotic_string(writer))\n",
              "├── write log: $(writer.write_log)\n",
              "└── Maximum file size: $(humanize_filesize(writer.max_filesize))"
              )
end

function save_diags_string(writer::PlanktonOutputWriter)
    if writer.save_diags
        return "diagnostics are saved as $(writer.diags_file)"
    else
        return "diagnostics are not saved"
    end
end
function save_phyto_string(writer::PlanktonOutputWriter)
    if writer.save_phytoplankton
        return string("phytoplankton is saved as $(writer.phytoplankton_file)\n",
                      "│   ├── every $(writer.phytoplankton_iteration_interval) seconds\n",
                      "│   └── including: $(writer.phytoplankton_include)")
    else
        return "phytoplankton are not saved"
    end
end
function save_abiotic_string(writer::PlanktonOutputWriter)
    if writer.save_abiotic_particle
        return string("abiotic particle is saved as $(writer.abiotic_particle_file)\n",
                      "│   ├── every $(writer.abiotic_particle_iteration_interval) seconds\n",
                      "│   └── including: $(writer.abiotic_particle_include)")
    else
        return "abiotic particle are not saved"
    end
end
function humanize_filesize(s::Number)
    suffix = ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"]
    biggest_suffix = suffix[1]
    value = s
    unit = 1024.0f0
    for power in 1:length(suffix)
        unit = 1024^power
        biggest_suffix = suffix[power]
        value < unit && break
    end
    s = 1024 * value / unit
    return @sprintf("%3.1f %s", s, biggest_suffix)
end
