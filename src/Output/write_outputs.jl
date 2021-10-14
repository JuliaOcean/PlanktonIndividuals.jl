##### write model outputs
function write_output!(writer::PlanktonOutputWriter, model::PlanktonModel, diags::PlanktonDiagnostics, ΔT)
    if isa(writer, Nothing)
        return nothing
    else
        if writer.write_log
            write_species_dynamics(model.t, model.individuals.phytos,
                                   writer.filepath, model.mode)
        end

        if writer.save_diags
            if model.t % (diags.time_interval) == 0.0
                if filesize(writer.diags_file) ≥ writer.max_filesize
                    start_next_diags_file(writer)
                end
                write_diags_to_jld2(diags, writer.diags_file, model.t,
                                    diags.time_interval/ΔT, model.grid)
            end
        end

        if writer.save_plankton
                if filesize(writer.plankton_file) ≥ writer.max_filesize
                    start_next_plankton_file(writer)
                end
            write_individuals_to_jld2(model.individuals.phytos, writer.plankton_file, model.t,
                                      writer.plankton_include)
        end
    end
end

function start_next_diags_file(writer::PlanktonOutputWriter)
    if writer.part_diags == 1
        part1_path = replace(writer.diags_file, r".jld2$" => "_part1.jld2")
        mv(writer.diags_file, part1_path, force=writer.force)
        writer.diags_file = part1_path
    end

    writer.part_diags += 1
    writer.diags_file = replace(writer.diags_file, r"part\d+.jld2$" => "part" * string(writer.part_diags) * ".jld2")
end

function start_next_plankton_file(writer::PlanktonOutputWriter)
    if writer.part_plankton == 1
        part1_path = replace(writer.plankton_file, r".jld2$" => "_part1.jld2")
        mv(writer.plankton_file, part1_path, force=writer.force)
        writer.plankton_file = part1_path
    end

    writer.part_plankton += 1
    writer.plankton_file = replace(writer.plankton_file, r"part\d+.jld2$" => "part" * string(writer.part_plankton) * ".jld2")
end

##### write a brief summary of each species at each time step into a txt file
function write_species_dynamics(t::Int64, phytos, filepath, mode::QuotaMode)
    for i in 1:length(phytos)
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        size_ave=  dot(phytos[i].data.Sz,  phytos[i].data.ac) / pop
        Bm_ave  =  dot(phytos[i].data.Bm,  phytos[i].data.ac) / pop
        Cq_ave  =  dot(phytos[i].data.Cq,  phytos[i].data.ac) / pop
        Nq_ave  =  dot(phytos[i].data.Nq,  phytos[i].data.ac) / pop
        Pq_ave  =  dot(phytos[i].data.Pq,  phytos[i].data.ac) / pop
        Chl_ave =  dot(phytos[i].data.chl, phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(filepath*"dynamic_species"*lpad(i,3,"0")*".txt","a");
        println(io,@sprintf("%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E",
                            day,hour,pop,gen_ave,age_ave,size_ave,Bm_ave,Cq_ave,Nq_ave,Pq_ave,Chl_ave))
        close(io);
    end
end
function write_species_dynamics(t::Int64, phytos, filepath, mode::CarbonMode)
    for i in 1:length(phytos)
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        size_ave=  dot(phytos[i].data.Sz,  phytos[i].data.ac) / pop
        Bm_ave  =  dot(phytos[i].data.Bm,  phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(filepath*"dynamic_species"*lpad(i,3,"0")*".txt","a");
        println(io,@sprintf("%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E",
                            day,hour,pop,gen_ave,age_ave,size_ave,Bm_ave))
        close(io);
    end
end

"""
    write_individuals_to_bin(phytos, filepath, t)
write model output of individuals at each time step to a binary file

Keyword Arguments
=================
- `phytos`: `NamedTuple` of a list of `individual` species.
- `filepath`: The file path to store JLD2 files.
- `t`: Current time of `model` in second, usually starting from 0.
- `atts` (optional): attributes of individuals to save, default `(:x, :y, :z)`
"""
function write_individuals_to_jld2(phytos::NamedTuple, filepath, t, atts)
    jldopen(filepath, "a+") do file
        for sp in keys(phytos)
            spi = NamedTuple{atts}([getproperty(phytos[sp].data, att) for att in atts])
            for att in atts
                file[lpad(t, 10, "0")*"/"*string(sp)*"/"*string(att)] = Array(spi[att])
            end
        end
    end
end

"""
    write_diags_to_jld2(diags, filepath, t, ncounts, grid)
write model output of individuals at each time step to a binary file

Keyword Arguments
=================
- `diags`: `NamedTuple` of a list of diagnostics at current time step.
- `filepath`: The file path to store JLD2 files.
- `t`: Current time of `model` in second, usually starting from 0.
- `ncounts`: the number of time steps included in each diagnostic
- `grid`: grid information used to exclude halo points.
"""
function write_diags_to_jld2(diags, filepath, t, ncounts, grid)
    jldopen(filepath, "a+") do file
        for key in keys(diags.tracer)
            file[lpad(t, 10, "0")*"/nut/"*string(key)] =
                Array(interior(diags.tracer[key], grid)) ./ ncounts
        end
        for sp in keys(diags.plankton)
            for proc in keys(diags.plankton[sp])
                file[lpad(t, 10, "0")*"/"*string(sp)*"/"*string(proc)] =
                    Array(interior(diags.plankton[sp][proc],grid)) ./ ncounts 
            end
        end
    end
    ##### zeros diags
    for tr in diags.tracer
        tr .= 0.0
    end
    for sp in diags.plankton
        for proc in sp
            proc .= 0.0
        end
    end
end