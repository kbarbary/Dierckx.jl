module Dierckx

export GridSpline,
       evaluate,
       evalgrid

const ddierckx = joinpath(dirname(@__FILE__),
                          "../deps/src/ddierckx/libddierckx")

# Ensure library is available.
if (dlopen_e(ddierckx) == C_NULL)
    error("Dierckx not properly installed. Run Pkg.build(\"Dierckx\")")
end

const _fit_messages = [1=>"""
The required storage space exceeds the available storage space: nxest
or nyest too small, or s too small.
The weighted least-squares spline corresponds to the current set of
knots.""",
                       2=>"""
A theoretically impossible result was found during the iteration
process for finding a smoothing spline with fp = s: s too small or
badly chosen eps.
Weighted sum of squared residuals does not satisfy abs(fp-s)/s < tol.""",
3=>"""
the maximal number of iterations maxit (set to 20 by the program)
allowed for finding a smoothing spline with fp=s has been reached:
s too small.
Weighted sum of squared residuals does not satisfy abs(fp-s)/s < tol.""",
4=>"""
No more knots can be added because the number of b-spline coefficients
(nx-kx-1)*(ny-ky-1) already exceeds the number of data points m:
either s or m too small.
The weighted least-squares spline corresponds to the current set of
knots.""",
5=>"""
No more knots can be added because the additional knot would (quasi)
coincide with an old one: s too small or too large a weight to an
inaccurate data point.
The weighted least-squares spline corresponds to the current set of
knots.""",
10=>"""
Error on entry, no approximation returned. The following conditions
must hold:
xb<=x[i]<=xe, yb<=y[i]<=ye, w[i]>0, i=0..m-1
If iopt==-1, then
  xb<tx[kx+1]<tx[kx+2]<...<tx[nx-kx-2]<xe
  yb<ty[ky+1]<ty[ky+2]<...<ty[ny-ky-2]<ye""",
  -3=>"""
The coefficients of the spline returned have been computed as the
minimal norm least-squares solution of a (numerically) rank deficient
system (deficiency=%i). If deficiency is large, the results may be
inaccurate. Deficiency may strongly depend on the value of eps."""
]

const _evaluate_message = """Invalid input data. Restrictions:
length(x) >= 1, length(y) >= 1
x(i-1) <= x(i), i=2,...,length(x)
y(j-1) <= y(j), j=2,...,length(y)
"""

# Currently only Float64, will add Float32 when we add
# single-precision dierckx library.
type GridSpline{T<:Union(Float64)}
    tx::Vector{T}
    ty::Vector{T}
    c::Vector{T}
    fp::T
    kx::Int
    ky::Int
end

# size of z should be (nx, ny).
function GridSpline(x::Vector{Float64}, y::Vector{Float64},
                    z::Array{Float64,2}; kx::Int=3, ky::Int=3,
                    s::Float64=0.0)
    mx = length(x)
    my = length(y)
    assert(size(z, 1) == mx && size(z, 2) == my)

    # Bounds
    xb = x[1]
    xe = x[end]
    yb = y[1]
    ye = y[end]
    nxest = mx+kx+1
    nyest = my+ky+1

    assert(mx > kx)
    assert(my > ky)
    assert(nxest >= 2 * (kx + 1))
    assert(nyest >= 2 * (ky + 1))

    # Return values
    nx = Array(Int32, 1)
    tx = Array(Float64, nxest)
    ny = Array(Int32, 1)
    ty = Array(Float64, nyest)
    c = Array(Float64, (nxest-kx-1) * (nyest-ky-1))
    fp = Array(Float64, 1)

    # work arrays
    lwrk = (4 + nyest * (mx+2*ky+5) + nxest * (2*kx+5) +
            my*(ky+1) + mx*(kx+1) + max(mx, nyest))
    wrk = Array(Float64, lwrk)
    kwrk = 3 + mx + my + nxest + nyest
    iwrk = Array(Int32, kwrk)

    ier = Array(Int32, 1)

    # Note regarding argument order:
    # The fortran function expects z to have shape (my, mx)
    # but we'd rather have x be the fast axis in z.
    # So, all the x and y related inputs are switched here from
    # what the fortran documentation says.
    ccall((:regrid_, ddierckx), Void,
          (Ptr{Int32},  # iopt
           Ptr{Int32}, Ptr{Float64},  # my, y
           Ptr{Int32}, Ptr{Float64},  # mx, x
           Ptr{Float64},  # z
           Ptr{Float64}, Ptr{Float64},  # yb, ye
           Ptr{Float64}, Ptr{Float64},  # xb, xe
           Ptr{Int32}, Ptr{Int32}, Ptr{Float64},  # ky, kx, s
           Ptr{Int32}, Ptr{Int32},  # nyest, nxest
           Ptr{Int32}, Ptr{Float64},  # ny, ty
           Ptr{Int32}, Ptr{Float64},  # nx, tx
           Ptr{Float64}, Ptr{Float64},  # c, fp
           Ptr{Float64}, Ptr{Int32},  # wrk, lwrk
           Ptr{Int32}, Ptr{Int32},  # iwrk, lwrk
           Ptr{Int32}),  # ier
          &0f0, &my, y, &mx, x, z, &yb, &ye, &xb, &xe, &ky, &kx, &s,
          &nyest, &nxest, ny, ty, nx, tx, c, fp, wrk, &lwrk, iwrk, &kwrk,
          ier)

    # Resize output arrays to the size actually used.
    resize!(tx, nx[1])
    resize!(ty, ny[1])
    resize!(c, (nx[1] - kx - 1) * (ny[1] - ky - 1))

    if !(ier[1] == 0 || ier[1] == -1 || ier[1] == -2)
        msg = _fit_messages[ier[1]]
        error(msg)
    end

    return GridSpline(tx, ty, c, fp[1], kx, ky)
end


function evaluate(spline::GridSpline, x::Vector{Float64},
                  y::Vector{Float64})
    m = length(x)
    assert(length(y) == m)

    ier = Array(Int32, 1)  # errors
    lwrk = spline.kx + spline.ky + 2
    wrk = Array(Float64, lwrk)
    z = Array(Float64, m)  # Return values

    ccall((:bispeu_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # ty, ny
           Ptr{Float64}, Ptr{Int32},  # tx, nx
           Ptr{Float64},  # c
           Ptr{Int32}, Ptr{Int32},  # ky, kx
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # y, x, z
           Ptr{Int32},  # m
           Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),  # wrk, lwrk, ier
          spline.ty, &length(spline.ty), spline.tx, &length(spline.tx),
          spline.c, &spline.ky, &spline.kx, y, x, z, &m, wrk, &lwrk, ier)

    if ier[1] != 0
        error(_evaluate_message)
    end

    return z
end

# Evaluate spline on the grid spanned by the input arrays.
function evalgrid(spline::GridSpline, x::Vector{Float64}, y::Vector{Float64})

    mx = length(x)
    my = length(y)
    ier = Array(Int32, 1)
    lwrk = mx*(spline.kx + 1) + my*(spline.ky + 1)
    wrk = Array(Float64, lwrk)
    kwrk = mx * my
    iwrk = Array(Int32, kwrk)

    z = Array(Float64, mx, my)  # Return values

    ccall((:bispev_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # ty, ny
           Ptr{Float64}, Ptr{Int32},  # tx, nx
           Ptr{Float64},  # c
           Ptr{Int32}, Ptr{Int32},  # ky, kx
           Ptr{Float64}, Ptr{Int32},  # y, my
           Ptr{Float64}, Ptr{Int32},  # x, mx
           Ptr{Float64},  # z
           Ptr{Float64}, Ptr{Int32},  # wrk, lwrk
           Ptr{Int32}, Ptr{Int32},  # iwrk, kwrk
           Ptr{Int32}),  # ier
          spline.ty, &length(spline.ty), spline.tx, &length(spline.tx),
          spline.c, &spline.ky, &spline.kx, y, &my, x, &mx, z,
          wrk, &lwrk, iwrk, &kwrk, ier)

    if ier[1] != 0
        error(_evalute_message)
    end

    return z
end


end # module
