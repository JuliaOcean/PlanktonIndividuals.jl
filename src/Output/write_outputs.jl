##### write model outputs
function write_output!(writer::Union{PlanktonOutputWriter, Nothing}, model::PlanktonModel, diags::Union{PlanktonDiagnostics,Nothing}, ΔT)
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
                write_diags_to_jld2(diags, writer.diags_file, model.t, model.iteration,
                                    diags.time_interval/ΔT, model.grid)
            end
        end

        if writer.save_plankton
            if model.t % (writer.plankton_time_interval) == 0.0
                if filesize(writer.plankton_file) ≥ writer.max_filesize
                    start_next_plankton_file(writer)
                end
                write_individuals_to_jld2(model.individuals.phytos, writer.plankton_file, model.t,
                                          model.iteration, writer.plankton_include)
            end
        end
    end
end

function start_next_diags_file(writer::PlanktonOutputWriter)
    if writer.part_diags == 1
        part1_path = replace(writer.diags_file, r".jld2$" => "_part1.jld2")
        mv(writer.diags_file, part1_path, force=true)
        writer.diags_file = part1_path
    end

    writer.part_diags += 1
    writer.diags_file = replace(writer.diags_file, r"part\d+.jld2$" => "part" * string(writer.part_diags) * ".jld2")
end

function start_next_plankton_file(writer::PlanktonOutputWriter)
    if writer.part_plankton == 1
        part1_path = replace(writer.plankton_file, r".jld2$" => "_part1.jld2")
        mv(writer.plankton_file, part1_path, force=true)
        writer.plankton_file = part1_path
    end

    writer.part_plankton += 1
    writer.plankton_file = replace(writer.plankton_file, r"part\d+.jld2$" => "part" * string(writer.part_plankton) * ".jld2")
end

##### write a brief summary of each species at each time step into a txt file
function write_species_dynamics(t::Int64, phytos, filepath, mode::MacroMolecularMode)
    for i in 1:length(phytos)
        file = joinpath(filepath, "dynamic_species"*lpad(i,3,"0")*".txt")
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        PRO_ave =  dot(phytos[i].data.PRO, phytos[i].data.ac) / pop
        DNA_ave =  dot(phytos[i].data.DNA, phytos[i].data.ac) / pop
        RNA_ave =  dot(phytos[i].data.RNA, phytos[i].data.ac) / pop
        CH_ave  =  dot(phytos[i].data.CH,  phytos[i].data.ac) / pop
        NST_ave =  dot(phytos[i].data.NST, phytos[i].data.ac) / pop
        PST_ave =  dot(phytos[i].data.PST, phytos[i].data.ac) / pop
        Chl_ave =  dot(phytos[i].data.Chl, phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(file,"a");
        println(io,@sprintf("%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E  %.8E  %.8E",
                            day,hour,pop,gen_ave,age_ave,PRO_ave,DNA_ave,RNA_ave,CH_ave,NST_ave,PST_ave,Chl_ave))
        close(io);
    end
end
function write_species_dynamics(t::Int64, phytos, filepath, mode::QuotaMode)
    for i in 1:length(phytos)
        file = joinpath(filepath, "dynamic_species"*lpad(i,3,"0")*".txt")
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        size_ave=  dot(phytos[i].data.Sz,  phytos[i].data.ac) / pop
        Bm_ave  =  dot(phytos[i].data.Bm,  phytos[i].data.ac) / pop
        Cq_ave  =  dot(phytos[i].data.Cq,  phytos[i].data.ac) / pop
        Nq_ave  =  dot(phytos[i].data.Nq,  phytos[i].data.ac) / pop
        Pq_ave  =  dot(phytos[i].data.Pq,  phytos[i].data.ac) / pop
        Chl_ave =  dot(phytos[i].data.Chl, phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(file,"a");
        println(io,@sprintf("%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E",
                            day,hour,pop,gen_ave,age_ave,size_ave,Bm_ave,Cq_ave,Nq_ave,Pq_ave,Chl_ave))
        close(io);
    end
end
function write_species_dynamics(t::Int64, phytos, filepath, mode::CarbonMode)
    for i in 1:length(phytos)
        file = joinpath(filepath, "dynamic_species"*lpad(i,3,"0")*".txt")
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        size_ave=  dot(phytos[i].data.Sz,  phytos[i].data.ac) / pop
        Bm_ave  =  dot(phytos[i].data.Bm,  phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(file,"a");
        println(io,@sprintf("%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E",
                            day,hour,pop,gen_ave,age_ave,size_ave,Bm_ave))
        close(io);
    end
end

function write_individuals_to_jld2(phytos::NamedTuple, filepath, t, iter, atts)
    jldopen(filepath, "a+") do file
        file["timeseries/t/$iter"] = t
        for sp in keys(phytos)
            spi = NamedTuple{atts}([getproperty(phytos[sp].data, att) for att in atts])
            for att in atts
                file["timeseries/$sp/$att/$iter"] = Array(spi[att])
            end
        end
    end
end

function write_diags_to_jld2(diags, filepath, t, iter, ncounts, grid)
    jldopen(filepath, "a+") do file
        file["timeseries/t/$iter"] = t
        for key in keys(diags.tracer)
            file["timeseries/$key/$iter"] = Array(interior(diags.tracer[key], grid)) ./ ncounts
        end
        for sp in keys(diags.plankton)
            for proc in keys(diags.plankton[sp])
                file["timeseries/$sp/$proc/$iter"] = Array(interior(diags.plankton[sp][proc],grid)) ./ ncounts 
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
