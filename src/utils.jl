##### copy external velocities into the model
function vel_copy!(vel::NamedTuple, u, v, w, arch::Architecture, g::Grids)
    vel.u.data[g.x⁻:g.x⁺, g.y⁻:g.y⁺, g.z⁻:g.z⁺] .= u
    vel.v.data[g.x⁻:g.x⁺, g.y⁻:g.y⁺, g.z⁻:g.z⁺] .= v
    vel.w.data[g.x⁻:g.x⁺, g.y⁻:g.y⁺, g.z⁻:g.z⁺] .= w

    fill_halo_vel!(vel, g)
end

##### sum up nutrient consumption counts into nutrient tendencies
function cts_to_Gcs!(plk, cts, g::Grids)
    for i in 1:length(nut_names)
        @inbounds plk[i].data[g.x⁻:g.x⁺, g.y⁻:g.y⁺, g.z⁻:g.z⁺] .= sum(cts[:,:,:,:,i], dims=4)[:,:,:,1]
    end
end

"""
    sub_nut_tendency!(a, b, c)
subtract one tendency from total tendencies
"""
function sub_nut_tendency!(a, b, c)
    for i in 1:length(a)
        @inbounds a[i].data .= b[i].data .- c[i].data
    end
end

function zero_fields!(a)
    for i in 1:length(a)
        @inbounds a[i].data .= 0.0
    end
end
# """
#     add_nut_tendency!(a, b)
# add one tendency to another tendency
# """
# function add_nut_tendency!(a, b)
#     for i in 1:length(a)
#         a[i].data .= a[i].data .+ b[i].data
#     end
# end
