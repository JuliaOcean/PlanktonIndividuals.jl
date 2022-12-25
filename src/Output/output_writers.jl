mutable struct PlanktonOutputWriter
    filepath::String
    write_log::Bool
    save_diags::Bool
    save_plankton::Bool
    diags_file::String
    plankton_file::String
    plankton_include::Tuple
    plankton_iteration_interval::Int
    max_filesize::Number # in Bytes
    part_diags::Int
    part_plankton::Int
end

"""
    PlanktonOutputWriter(;dir = "./results",
                               diags_prefix = "diags",
                               plankton_prefix = "plankton",
                               write_log = false,
                               save_diags = false,
                               save_plankton = false,
                               plankton_include = (:x, :y, :z, :Sz),
                               plankton_iteration_interval = 1,
                               max_filesize = Inf,
                               )
Generate a `PlanktonOutputWriter` structure which includes settings for model outputs

Keyword Arguments (Optional)
============================
- `dir`: The directory to store model outputs, "./results" by default
- `diags_prefix`: Descriptive filename prefixed to diagnostic output files.
- `plankton_prefix`: Descriptive filename prefixed to plankton output files.
- `write_log`: write model logs which contain global averages of simulated plankton, default: `false`.
- `save_diags`: write diagnostics to disk, default: `false`.
- `save_plankton`: write plankton to disk, default: `false`.
- `plankton_include`: list of plankton properties to save, default: `(:x, :y, :z, :Sz)`.
- `plankton_iteration_interval`: The time interval that plankton are saved, 1 timestep by default.
- `max_filesize`: The writer will stop writing to the output file once the file size exceeds `max_filesize`,
                    and write to a new one with a consistent naming scheme ending in `part1`, `part2`, etc.
                    default: `Inf`.
"""
function PlanktonOutputWriter(;dir = "./results",
                               diags_prefix = "diags",
                               plankton_prefix = "plankton",
                               write_log = false,
                               save_diags = false,
                               save_plankton = false,
                               plankton_include = (:x, :y, :z, :Sz),
                               plankton_iteration_interval = 1,
                               max_filesize = Inf
                               )

    isdir(dir) && rm(dir, recursive=true)
    mkdir(dir)

    diags_file = ""
    plankton_file = ""

    if save_diags
        diags_file = joinpath(dir, diags_prefix*".jld2")
    end
    if save_plankton
        plankton_file = joinpath(dir, plankton_prefix*".jld2")
    end

    return PlanktonOutputWriter(dir, write_log, save_diags, save_plankton, diags_file, plankton_file,
                                plankton_include, plankton_iteration_interval, max_filesize, 1, 1)
end

function show(io::IO, writer::PlanktonOutputWriter)
    print(io, "PlanktonOutputWriter:\n",
              "├── files are saved at $(writer.filepath)\n",
              "├── $(save_diags_string(writer))\n",
              "├── $(save_inds_string(writer))\n",
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
function save_inds_string(writer::PlanktonOutputWriter)
    if writer.save_plankton
        return string("individuals are saved as $(writer.plankton_file)\n",
                      "│   ├── every $(writer.plankton_iteration_interval) seconds\n",
                      "│   └── including: $(writer.plankton_include)")
    else
        return "individuals are not saved"
    end
end
function humanize_filesize(s::Number)
    suffix = ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"]
    biggest_suffix = suffix[1]
    value = Float64(s)
    unit = 1024.0
    for power in 1:length(suffix)
        unit = 1024^power
        biggest_suffix = suffix[power]
        value < unit && break
    end
    s = 1024 * value / unit
    return @sprintf("%3.1f %s", s, biggest_suffix)
end
