"""
    default_PARF(grid, ΔT, iterations)
Generate default hourly surface PAR.
"""
function default_PARF(grid, ΔT, iterations)
    PAR = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3871666666666666, 87.10258333333333, 475.78150000000016, 929.2737916666669,
           1232.3633333333337, 1638.918916666667, 1823.7921666666664, 1906.2769583333336, 1776.0280416666667,
           1678.5026249999999, 1410.216666666667, 815.4129583333336, 525.104, 135.993, 2.9493750000000003, 0.0, 0.0, 0.0,]
    total_days = Int((ΔT * iterations) ÷ 86400 + 1)
    PAR_day = zeros(grid.Nx, grid.Ny, 24)
    for i in 1:24
        PAR_day[:,:,i] .= PAR[i]
    end
    PAR_domain = repeat(PAR_day, outer = (1,1,total_days))
    return PAR_domain
end

"""
    default_temperature(grid, ΔT, iterations)
Generate default hourly temperature.
"""
function default_temperature(grid, ΔT, iterations)
    temp = [26.646446609406727, 26.56698729810778, 26.517037086855467, 26.5, 26.517037086855467, 26.56698729810778,
            26.646446609406727, 26.75, 26.87059047744874, 27.0, 27.12940952255126, 27.25, 27.353553390593273,
            27.43301270189222, 27.482962913144533, 27.5, 27.482962913144533, 27.43301270189222, 27.353553390593273,
            27.25, 27.12940952255126, 27.0, 26.87059047744874, 26.75]
    total_days = Int((ΔT * iterations) ÷ 86400 + 1)
    temp_day = zeros(grid.Nx, grid.Ny, grid.Nz, 24)
    for i in 1:24
        temp_day[:,:,end,i] .= temp[i]
    end
    # vertical temperature gradient
    for j in grid.Nz-1:-1:1
        temp_day[:,:,j,:] .= temp_day[:,:,j+1,:] .+ 4.0e-4 * grid.zC[j]
    end
    temp_domain = repeat(temp_day, outer = (1,1,1,total_days))
    return temp_domain
end
