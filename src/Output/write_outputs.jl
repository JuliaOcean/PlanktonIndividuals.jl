##### write model outputs
function write_output!(writer::Union{PlanktonOutputWriter, Nothing}, model::PlanktonModel, diags::PlanktonDiagnostics)
    if isa(writer, Nothing)
        return nothing
    else
        if writer.write_log
            write_species_dynamics(model.t, model.individuals.phytos,
                                   writer.filepath, model.mode)
        end

        if writer.save_diags
            if model.iteration % diags.iteration_interval == 0.0f0
                if filesize(writer.diags_file) ≥ writer.max_filesize
                    start_next_diags_file(writer)
                end
                write_diags_to_jld2(diags, writer.diags_file, model.t, model.iteration,
                                    diags.iteration_interval, model.grid)
            end
        end

        if writer.save_phytoplankton
            if model.iteration % writer.phytoplankton_iteration_interval == 0.0f0
                if filesize(writer.phytoplankton_file) ≥ writer.max_filesize
                    start_next_phytoplankton_file(writer)
                end
                write_individuals_to_jld2(model.individuals.phytos, writer.phytoplankton_file, model.t,
                                          model.iteration, writer.phytoplankton_include)
            end
        end

        if writer.save_abiotic_particle
            if model.iteration % writer.abiotic_particle_iteration_interval == 0.0f0
                if filesize(writer.abiotic_particle_file) ≥ writer.max_filesize
                    start_next_abiotic_particle_file(writer)
                end
                write_individuals_to_jld2(model.individuals.abiotics, writer.abiotic_particle_file, model.t,
                                          model.iteration, writer.abiotic_particle_include)
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

function start_next_phytoplankton_file(writer::PlanktonOutputWriter)
    if writer.part_phytoplankton == 1
        part1_path = replace(writer.phytoplankton_file, r".jld2$" => "_part1.jld2")
        mv(writer.phytoplankton_file, part1_path, force=true)
        writer.phytoplankton_file = part1_path
    end

    writer.part_phytoplankton += 1
    writer.phytoplankton_file = replace(writer.phytoplankton_file, r"part\d+.jld2$" => "part" * string(writer.part_phytoplankton) * ".jld2")
end

function start_next_abiotic_particle_file(writer::PlanktonOutputWriter)
    if writer.part_abiotic_particle == 1
        part1_path = replace(writer.abiotic_particle_file, r".jld2$" => "_part1.jld2")
        mv(writer.abiotic_particle_file, part1_path, force=true)
        writer.abiotic_particle_file = part1_path
    end

    writer.part_abiotic_particle += 1
    writer.abiotic_particle_file = replace(writer.abiotic_particle_file, r"part\d+.jld2$" => "part" * string(writer.part_abiotic_particle) * ".jld2")
end

##### write a brief summary of each species at each time step into a txt file
function write_species_dynamics(t::AbstractFloat, phytos, filepath, mode::MacroMolecularMode)
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
function write_species_dynamics(t::AbstractFloat, phytos, filepath, mode::QuotaMode)
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
function write_species_dynamics(t::AbstractFloat, phytos, filepath, mode::IronEnergyMode)
    for i in 1:length(phytos)
        file = joinpath(filepath, "dynamic_species"*lpad(i,3,"0")*".txt")
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        size_ave=  dot(phytos[i].data.Sz,  phytos[i].data.ac) / pop
        Bm_ave  =  dot(phytos[i].data.Bm,  phytos[i].data.ac) / pop
        En_ave  =  dot(phytos[i].data.En,  phytos[i].data.ac) / pop
        CH_ave  =  dot(phytos[i].data.CH,  phytos[i].data.ac) / pop
        qNO3_ave=  dot(phytos[i].data.qNO3,  phytos[i].data.ac) / pop
        qNH4_ave=  dot(phytos[i].data.qNH4,  phytos[i].data.ac) / pop
        qP_ave  =  dot(phytos[i].data.qP,  phytos[i].data.ac) / pop
        qFe_ave =  dot(phytos[i].data.qFe,  phytos[i].data.ac) / pop
        Chl_ave =  dot(phytos[i].data.Chl, phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(file,"a");
        println(io,@sprintf(
            "%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E  %.8E  %.8E  %.8E",
            day,hour,pop,gen_ave,age_ave,size_ave,Bm_ave,En_ave,CH_ave,qNO3_ave,qNH4_ave,qP_ave,qFe_ave,Chl_ave))
        close(io);
    end
end
function write_species_dynamics(t::AbstractFloat, phytos, filepath, mode::CarbonMode)
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

function write_individuals_to_jld2(particles::NamedTuple, filepath, t, iter, atts)
    jldopen(filepath, "a+") do file
        file["timeseries/t/$iter"] = t
        for sp in keys(particles)
            spi = NamedTuple{atts}([getproperty(particles[sp].data, att) for att in atts])
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
        for sp in keys(diags.phytos)
            for proc in keys(diags.phytos[sp])
                file["timeseries/phyto/$sp/$proc/$iter"] = Array(interior(diags.phytos[sp][proc],grid)) ./ ncounts 
            end
        end
        for sa in keys(diags.abiotics)
            for proc in keys(diags.abiotics[sa])
                file["timeseries/abiotic/$sa/$proc/$iter"] = Array(interior(diags.abiotics[sa][proc],grid)) ./ ncounts 
            end
        end
    end
    ##### zeros diags
    for tr in diags.tracer
        tr .= 0.0f0
    end
    for sp in diags.phytos
        for proc in sp
            proc .= 0.0f0
        end
    end
    for sp in diags.abiotics
        for proc in sp
            proc .= 0.0f0
        end
    end
end
