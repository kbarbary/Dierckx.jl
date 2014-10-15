module Dierckx

const libpath = joinpath(dirname(@__FILE__),
                         "../deps/src/ddierckx/libddierckx")
ddierckx = dlopen(libpath)

function regrid(x, y, z)

    mx::Int32 = length(x)
    my::Int32 = length(y)
    assert(size(z, 1) == mx && size(z, 2) == my)

    xb = x[1]
    xe = x[end]
    yb = y[1]
    ye = y[end]
    kx::Int32 = 3
    ky::Int32 = 3
    s = 0.0
    nxest::Int32 = mx+kx+1
    nyest::Int32 = my+ky+1

    assert(mx > kx)
    assert(my > ky)
    assert(nxest >= 2 * (kx + 1))
    assert(nxest >= 2 * (kx + 1))

    nx = Array(Int32, 1)
    tx = Array(Float64, nxest)
    ny = Array(Int32, 1)
    ty = Array(Float64, nyest)
    c = Array(Float64, (nxest-kx-1) * (nyest-ky-1))
    fp = Array(Float64, 1)

    # size of work array
    lwrk::Int32 = (4 + nxest * (my+2*kx+5) + nyest * (2*ky+5) +
            mx*(kx+1) + my*(ky+1) + max(my, nxest))
    wrk = Array(Float64, lwrk)
    kwrk::Int32 = 3 + mx + my + nxest + nyest
    iwrk = Array(Int32, kwrk)

    ier = Array(Int32, 1)

    ccall((:regrid_, libpath), Void,
          (Ptr{Int32},  # iopt
           Ptr{Int32}, Ptr{Float64},  # mx, x
           Ptr{Int32}, Ptr{Float64},  # my, y
           Ptr{Float64},  # z
           Ptr{Float64}, Ptr{Float64},  # xb, xe
           Ptr{Float64}, Ptr{Float64},  # yb, ye
           Ptr{Int32}, Ptr{Int32}, Ptr{Float64},  # kx, ky, s
           Ptr{Int32}, Ptr{Int32},  # nxest, nyest
           Ptr{Int32}, Ptr{Float64},  # nx, tx
           Ptr{Int32}, Ptr{Float64},  # ny, ty
           Ptr{Float64}, Ptr{Float64},  # c, fp
           Ptr{Float64}, Ptr{Int32},  # wrk, lwrk
           Ptr{Int32}, Ptr{Int32},  # iwrk, lwrk
           Ptr{Int32}),  # ier
          &0f0, &mx, x, &my, y, z, &xb, &xe, &yb, &ye, &kx, &ky, &s,
          &nxest, &nyest, nx, tx, ny, ty, c, fp,
          wrk, &lwrk, iwrk, &kwrk, ier)

    # Resize output arrays to the size actually used.
    resize!(tx, nx[1])
    resize!(ty, ny[1])
    resize!(c, (nx[1] - kx - 1) * (ny[1] - ky - 1))

    if ier[1] >= 1
        error("regrid error $(ier[1])")
    end

    return tx, ty, c, fp
end

end # module
