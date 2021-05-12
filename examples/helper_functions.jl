
"""
    (uvels, vvels, wvels) = streamfunction_xz()

Define a two-dimensional flow field in the x-z directions from a simple streamfunction.
"""
function streamfunction_xz()
    scal = 2e-1
    f(x, y, z) = scal*sin(x*2π/128)*sin(z*π/128) #stream function

    ϕcorners=[f(x,0.,z) for x in 0:128, z in -128:0]
    ϕcenters=[f(x,0.,z) for x in 0.5:128, z in -128:-0.5]

    uu=-diff(ϕcorners,dims=2)[1:end-1,:]
    ww=diff(ϕcorners,dims=1)
    uu=reshape(uu,(128,1,128))
    ww=reshape(ww,(128,1,129))
    uu = reverse(uu, dims=3)
    ww = reverse(ww, dims=3)

    uvels = fill(uu, 2)
    vvels = fill(0*uu, 2)
    wvels = fill(ww, 2)
    uvels = cat(uvels..., dims=4)
    vvels = cat(vvels..., dims=4)
    wvels = cat(wvels..., dims=4)

    return uvels, vvels, wvels, ϕcenters
end
