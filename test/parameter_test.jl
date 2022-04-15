using PlanktonIndividuals.Parameters

function test_parameters()
    tmp_bgc = Dict("κh" => 1e-5, "κv" => 1e-5)
    param_bgc = update_bgc_params(tmp_bgc)
    @test param_bgc["κh"] == 1e-5
    @test param_bgc["κv"] == 1e-5

    tmp_phyt_carbon = Dict("PCmax" => [1e-5])
    param_phyt_carbon = update_phyt_params(tmp_phyt_carbon; N = 1, mode = CarbonMode())
    @test param_phyt_carbon["PCmax"] == [1e-5]

    tmp_phyt_quota = Dict("PCmax" => [1e-5, 1e-5, 1e-5])
    param_phyt_quota = update_phyt_params(tmp_phyt_quota; N = 3, mode = QuotaMode())
    @test param_phyt_quota["PCmax"] == [1e-5, 1e-5, 1e-5]

    return nothing
end

@testset "Parameterss" begin
    test_parameters()
end
