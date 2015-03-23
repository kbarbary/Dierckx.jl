module Dierckx

using Compat

export Spline1D,
       Spline2D,
       evaluate,
       derivative,
       evalgrid,
       get_knots,
       get_coeffs,
       get_residual

import Base: show

unixpath = "../deps/src/ddierckx/libddierckx"
winpath = "../deps/bin$WORD_SIZE/libddierckx"
const ddierckx = joinpath(dirname(@__FILE__), @unix? unixpath : winpath)

function __init__()
    # Ensure library is available.
    if (Libdl.dlopen_e(ddierckx) == C_NULL)
        error("Dierckx not properly installed. Run Pkg.build(\"Dierckx\")")
    end
end

# ----------------------------------------------------------------------------
# 1-d splines

const _fit1d_messages = @compat Dict(
2=>
"""A theoretically impossible result was found during the iteration
process for finding a smoothing spline with fp = s: s too small.
There is an approximation returned but the corresponding weighted sum
of squared residuals does not satisfy the condition abs(fp-s)/s <
tol.""",
3=>
"""The maximal number of iterations maxit (set to 20 by the program)
allowed for finding a smoothing spline with fp=s has been reached: s
too small. There is an approximation returned but the corresponding
weighted sum of squared residuals does not satisfy the condition
abs(fp-s)/s < tol.""",
10=>
"""Error on entry, no approximation returned. The following conditions
must hold:
1<=k<=5
x[1] < x[2] < ... < x[end]
w[i] > 0.0 for all i

Additionally, if spline knots are given:
length(xknots) <= length(x) + k + 1
x[1] < xknots[1] < xknots[k+2] < ... < xknots[end] < x[end]
The schoenberg-whitney conditions: there must be a subset of data points
xx[j] such that t[j] < xx[j] < t[j+k+1] for j=1,2,...,n-k-1""")


const _eval1d_messages = @compat Dict(
1=>
"""Input point out of range""",
10=>
"""Invalid input data. The following conditions must hold:
length(x) != 0 and xb <= x[1] <= x[2] <= ... x[end] <= xe""")

_translate_bc(bc::String) = (bc == "extrapolate" ? 0 :
                             bc == "zero" ? 1 :
                             bc == "error" ? 2 :
                             bc == "nearest" ? 3 :
                             error("unknown boundary condition: \"$bc\""))
_translate_bc(bc::Int) = (bc == 0 ? "extrapolate" :
                          bc == 1 ? "zero" :
                          bc == 2 ? "error" :
                          bc == 3 ? "nearest" : "")

type Spline1D
    t::Vector{Float64}
    c::Vector{Float64}
    k::Int
    bc::Int
    fp::Float64
end

get_knots(spl::Spline1D) = spl.t[spl.k+1:end-spl.k]
get_coeffs(spl::Spline1D) = spl.c[1:end-spl.k+1]
get_residual(spl::Spline1D) = spl.fp

function reallycompact(a::Vector)
    io = IOBuffer()
    if length(a) <= 5
        show(io, a)
    else
        write(io, "[")
        showcompact(io, a[1])
        write(io, ",")
        showcompact(io, a[2])
        write(io, " \u2026 ")
        showcompact(io, a[end-1])
        write(io, ",")
        showcompact(io, a[end])
        write(io, "]")
        write(io, " ($(length(a)) elements)")
    end
    seekstart(io)
    readall(io)
end


function show(io::IO, spl::Spline1D)
    print(io, """Spline1D(knots=$(reallycompact(get_knots(spl))), k=$(spl.k), extrapolation=\"$(_translate_bc(spl.bc))\", residual=$(spl.fp))""")
end

function Spline1D(x::Vector{Float64}, y::Vector{Float64};
                  w::Vector{Float64}=ones(length(x)),
                  k::Int=3, s::Real=0.0, bc::String="nearest")
    m = length(x)
    length(y) == m || error("length of x and y must match")
    length(w) == m || error("length of x and w must match")
    m > k || error("length(x) must be greater than k")

    # outputs
    nest = m+k+1
    n = Array(Int32, 1)
    t = Array(Float64, nest)
    c = Array(Float64, nest)
    fp = Array(Float64, 1)
    ier = Array(Int32, 1)

    # working space
    lwrk = m * (k+1) + nest * (7 + 3*k)
    wrk = Array(Float64, lwrk)
    iwrk = Array(Int32, nest)

    ccall((:curfit_, ddierckx), Void,
          (Ptr{Int32}, Ptr{Int32},  # iopt, m
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # x, y, w
           Ptr{Float64}, Ptr{Float64},  # xb, xe
           Ptr{Int32}, Ptr{Float64},  # k, s
           Ptr{Int32}, Ptr{Int32}, # nest, n
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # t, c, fp
           Ptr{Float64}, Ptr{Int32}, Ptr{Int32},  # wrk, lwrk, iwrk
           Ptr{Int32}),  # ier
          &0, &m, x, y, w, &x[1], &x[end], &k, &@compat(Float64(s)), &nest,
          n, t, c, fp, wrk, &lwrk, iwrk, ier)

    ier[1] <= 0 || error(_fit1d_messages[ier[1]])

    # resize output arrays
    resize!(t, n[1])
    resize!(c, n[1] - k - 1)

    return Spline1D(t, c, k, _translate_bc(bc), fp[1])
end


# version with user-supplied knots
function Spline1D(x::Vector{Float64}, y::Vector{Float64},
                  xknots::Vector{Float64};
                  w::Vector{Float64}=ones(length(x)),
                  k::Int=3, bc::String="nearest")
    m = length(x)
    length(y) == m || error("length of x and y must match")
    length(w) == m || error("length of x and w must match")
    m > k || error("k must be less than length(x)")

    # x knots
    # (k+1) knots will be added on either end of interior knots.
    n = length(xknots) + 2 * (k+1)
    t = Array(Float64, n)  # All knots
    for i in 1:length(xknots)
        t[i+k+1] = xknots[i]
    end

    # outputs
    c = Array(Float64, n)
    fp = Array(Float64, 1)
    ier = Array(Int32, 1)

    # working space
    lwrk = m * (k+1) + n * (7 + 3*k)
    wrk = Array(Float64, lwrk)
    iwrk = Array(Int32, n)

    ccall((:curfit_, ddierckx), Void,
          (Ptr{Int32}, Ptr{Int32},  # iopt, m
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # x, y, w
           Ptr{Float64}, Ptr{Float64},  # xb, xe
           Ptr{Int32}, Ptr{Float64},  # k, s
           Ptr{Int32}, Ptr{Int32}, # nest, n
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # t, c, fp
           Ptr{Float64}, Ptr{Int32}, Ptr{Int32},  # wrk, lwrk, iwrk
           Ptr{Int32}),  # ier
          &(-1), &m, x, y, w, &x[1], &x[end], &k, &(-1.0), &n,
          &n, t, c, fp, wrk, &lwrk, iwrk, ier)

    ier[1] <= 0 || error(_fit1d_messages[ier[1]])
    resize!(c, n-k-1)

    return Spline1D(t, c, k, _translate_bc(bc), fp[1])
end


function evaluate(spline::Spline1D, x::Vector{Float64})
    m = length(x)
    y = Array(Float64, m)
    ier = Array(Int32, 1)
    ccall((:splev_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # t, n
           Ptr{Float64}, Ptr{Int32},  # c, k
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32},  # x, y, m
           Ptr{Int32}, Ptr{Int32}),  # e, ier
          spline.t, &length(spline.t), spline.c, &spline.k,
          x, y, &m, &spline.bc, ier)
    ier[1] == 0 || error(_eval1d_messages[ier[1]])
    return y
end


function evaluate(spline::Spline1D, x::FloatingPoint)
    y = Array(Float64, 1)
    ier = Array(Int32, 1)
    ccall((:splev_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # t, n
           Ptr{Float64}, Ptr{Int32},  # c, k
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32},  # x, y, m
           Ptr{Int32}, Ptr{Int32}),  # e, ier
          spline.t, &length(spline.t), spline.c, &spline.k,
          &x, y, &1, &spline.bc, ier)
    ier[1] == 0 || error(_eval1d_messages[ier[1]])
    return y[1]
end


# TODO: should the function name be evalder, derivative, or grad?
#       or should it be integrated with evaluate, above?
#       (problem with that: derivative doesn't accept bc="nearest")
# TODO: should `nu` be `d`?
function derivative(spline::Spline1D, x::Vector{Float64}; nu::Integer=1)

    # Unlike evaluate, this Fortran function doesn't accept bc="nearest".
    # However, "nearest" gives a derivative of 0 for all nu outside the 
    # support range, so the answer will be the same as bc="zero".
    bc = spline.bc == 3 ? 1 : spline.bc

    # Check derivate requested.
    if !(1 <= nu <= spline.k)
        error("nu must be positive and less than or equal to spline order")
    end

    m = length(x)
    n = length(spline.t)
    wrk = Array(Float64, n)
    y = Array(Float64, m)
    ier = Array(Int32, 1)
    ccall((:splder_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # t, n
           Ptr{Float64}, Ptr{Int32},  # c, k
           Ptr{Int32},  # nu
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32},  # x, y, m
           Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),  # e, wrk, ier
          spline.t, &n, spline.c, &spline.k, &nu,
          x, y, &m, &bc, wrk, ier)
    ier[1] == 0 || error(_eval1d_messages[ier[1]])
    return y
end

derivative(spline::Spline1D, x::FloatingPoint; nu::Integer=1) =
    derivative(spline, Float64[x]; nu=nu)[1]


# ----------------------------------------------------------------------------
# 2-d splines

# NOTE REGARDING ARGUMENT ORDER: In the "grid" version of the Spline2D
# constructor and evaluators, the fortran functions expects z to have
# shape (my, mx), but we'd rather have x be the fast axis in z.  So,
# in the ccall()s in these methods, all the x and y related inputs are
# swapped with regard to what the Fortran documentation says.

const _fit2d_messages = @compat Dict(
-3=>

"""The coefficients of the spline returned have been computed as the
minimal norm least-squares solution of a (numerically) rank deficient
system (deficiency=%i). If deficiency is large, the results may be
inaccurate. Deficiency may strongly depend on the value of eps.""",

1=>

"""The required storage space exceeds the available storage space:
nxest or nyest too small, or s too small. Try increasing s.""",
# The weighted least-squares spline corresponds to the current set of knots.

2=>

"""A theoretically impossible result was found during the iteration
process for finding a smoothing spline with fp = s: s too small or
badly chosen eps.  Weighted sum of squared residuals does not satisfy
abs(fp-s)/s < tol.""",

3=>

"""the maximal number of iterations maxit (set to 20 by the program)
allowed for finding a smoothing spline with fp=s has been reached: s
too small.  Weighted sum of squared residuals does not satisfy
abs(fp-s)/s < tol.""",

4=>

"""No more knots can be added because the number of b-spline
coefficients (nx-kx-1)*(ny-ky-1) already exceeds the number of data
points m: either s or m too small.  The weighted least-squares spline
corresponds to the current set of knots.""",

5=>

"""No more knots can be added because the additional knot would
(quasi) coincide with an old one: s too small or too large a weight to
an inaccurate data point.  The weighted least-squares spline
corresponds to the current set of knots.""",

10=>

"""Error on entry, no approximation returned. The following conditions
must hold:
xb<=x[i]<=xe, yb<=y[i]<=ye, w[i]>0, i=0..m-1
If iopt==-1, then
  xb<tx[kx+1]<tx[kx+2]<...<tx[nx-kx-2]<xe
  yb<ty[ky+1]<ty[ky+2]<...<ty[ny-ky-2]<ye""")

const _eval2d_message = (
"""Invalid input data. Restrictions:
length(x) != 0, length(y) != 0
x[i-1] <= x[i] for i=2,...,length(x)
y[j-1] <= y[j] for j=2,...,length(y)
""")

type Spline2D
    tx::Vector{Float64}
    ty::Vector{Float64}
    c::Vector{Float64}
    kx::Int
    ky::Int
    fp::Float64
end

get_knots(spl::Spline2D) = (spl.tx[spl.kx+1:end-spl.kx],
                            spl.ty[spl.ky+1:end-spl.ky])
get_residual(spl::Spline2D) = spl.fp


# Helper function for calculating required size of work array in surfit.
# Note that x and y here are as in the Fortran documentation.
function calc_surfit_lwrk1(m, kx, ky, nxest, nyest)
    u = nxest - kx - 1
    v = nyest - ky - 1
    km = max(kx, ky) + 1
    ne = max(nxest, nyest)
    bx = kx*v + ky + 1
    by = ky*u + kx + 1
    b1 = b2 = 0
    if (bx<=by)
        b1 = bx
        b2 = bx + v - ky
    else
        b1 = by
        b2 = by + u - kx
    end
    return u*v*(2 + b1 + b2) + 2*(u+v+km*(m+ne)+ne-kx-ky) + b2 + 1
end

# Construct spline from unstructured data
function Spline2D(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64};
                  w::Vector{Float64}=ones(length(x)), kx::Int=3, ky::Int=3,
                  s::Real=0.0)

    # arary sizes
    m = length(x)
    (length(y) == length(z) == m) || error("lengths of x, y, z must match") 
    (length(w) == m) || error("length of w must match other inputs")

    nxest = max(kx+1+ceil(Int,sqrt(m/2)), 2*(kx+1))
    nyest = max(ky+1+ceil(Int,sqrt(m/2)), 2*(ky+1))
    nmax = max(nxest, nyest)

    eps = 1.0e-16

    # bounds
    xb = minimum(x)
    xe = maximum(x)
    yb = minimum(y)
    ye = maximum(y)

    # return values
    nx = Array(Int32, 1)
    tx = Array(Float64, nxest)
    ny = Array(Int32, 1)
    ty = Array(Float64, nyest)
    c = Array(Float64, (nxest-kx-1) * (nyest-ky-1))
    fp = Array(Float64, 1)
    ier = Array(Int32, 1)
    
    # work arrays
    # Note: in lwrk1, x and y are swapped on purpose.
    lwrk1 = calc_surfit_lwrk1(m, ky, kx, nyest, nxest)
    lwrk2 = 1
    kwrk = m + (nxest - 2*kx - 1) * (nyest - 2*ky - 1)
    wrk1 = Array(Float64, lwrk1)
    wrk2 = Array(Float64, lwrk2)
    iwrk = Array(Int32, kwrk)

    ccall((:surfit_, ddierckx), Void,
          (Ptr{Int32}, Ptr{Int32},  # iopt, m
           Ptr{Float64}, Ptr{Float64},  # y, x
           Ptr{Float64}, Ptr{Float64},  # z, w
           Ptr{Float64}, Ptr{Float64},  # yb, ye
           Ptr{Float64}, Ptr{Float64},  # xb, xe
           Ptr{Int32}, Ptr{Int32}, Ptr{Float64},  # ky, kx, s
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32},  # nyest, nxest, nmax
           Ptr{Float64},  # eps
           Ptr{Int32}, Ptr{Float64},  # ny, tx
           Ptr{Int32}, Ptr{Float64},  # nx, tx
           Ptr{Float64}, Ptr{Float64},  # c, fp
           Ptr{Float64}, Ptr{Int32},  # wrk1, lwrk1
           Ptr{Float64}, Ptr{Int32},  # wrk2, lwrk2
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),   # iwrk, kwrk, ier
          &0, &m, y, x, z, w, &yb, &ye, &xb, &xe, &ky, &kx, &@compat(Float64(s)),
          &nyest, &nxest, &nmax, &eps, ny, ty, nx, tx, c, fp,
          wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, ier)

    while ier[1] > 10
        # lwrk2 is too small, i.e., there is not enough workspace
        # for computing the minimal least-squares solution of a rank
        # deficient system of linear equations. ier gives the
        # requested value for lwrk2. Rerun with that value in "continue"
        # mode, with iopt = 1.
        lwrk2 = ier[1]
        resize!(wrk2, lwrk2)
        ccall((:surfit_, ddierckx), Void,
              (Ptr{Int32}, Ptr{Int32},  # iopt, m
               Ptr{Float64}, Ptr{Float64},  # y, x
               Ptr{Float64}, Ptr{Float64},  # z, w
               Ptr{Float64}, Ptr{Float64},  # yb, ye
               Ptr{Float64}, Ptr{Float64},  # xb, xe
               Ptr{Int32}, Ptr{Int32}, Ptr{Float64},  # ky, kx, s
               Ptr{Int32}, Ptr{Int32}, Ptr{Int32},  # nyest, nxest, nmax
               Ptr{Float64},  # eps
               Ptr{Int32}, Ptr{Float64},  # ny, tx
               Ptr{Int32}, Ptr{Float64},  # nx, tx
               Ptr{Float64}, Ptr{Float64},  # c, fp
               Ptr{Float64}, Ptr{Int32},  # wrk1, lwrk1
               Ptr{Float64}, Ptr{Int32},  # wrk2, lwrk2
               Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),   # iwrk, kwrk, ier
              &1, &m, y, x, z, w, &yb, &ye, &xb, &xe, &ky, &kx, &s,
              &nyest, &nxest, &nmax, &eps, ny, ty, nx, tx, c, fp,
              wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, ier)
    end

    if (ier[1] == 0 || ier[1] == -1 || ier[1] == -2)
        # good values, pass.
    elseif ier[1] < -2
        warn("""
        The coefficients of the spline returned have been
        computed as the minimal norm least-squares solution of a
        (numerically) rank deficient system. The rank is $(-ier[1]).
        The rank deficiency is $((nx[1]-kx-1)*(ny[1]-ky-1)+ier[1]).
        Especially if the rank deficiency is large the results may
        be inaccurate.""") 
        # "The results could also seriously depend on the value of
        # eps" (not in message because eps is currently not an input)
    else
        error(_fit2d_messages[ier[1]])
    end
            
    # Resize output arrays to the size actually used.
    resize!(tx, nx[1])
    resize!(ty, ny[1])
    resize!(c, (nx[1] - kx - 1) * (ny[1] - ky - 1))
    
    return Spline2D(tx, ty, c, kx, ky, fp[1])
end



# Construct spline from data on a grid.
function Spline2D(x::Vector{Float64}, y::Vector{Float64}, z::Array{Float64,2};
                  kx::Int=3, ky::Int=3, s::Real=0.0)
    mx = length(x)
    my = length(y)
    assert(size(z, 1) == mx && size(z, 2) == my)

    mx > kx || error("length(x) must be greater than kx")
    my > ky || error("length(y) must be greater than ky")
    
    # Bounds
    xb = x[1]
    xe = x[end]
    yb = y[1]
    ye = y[end]
    nxest = mx+kx+1
    nyest = my+ky+1

    # Return values
    nx = Array(Int32, 1)
    tx = Array(Float64, nxest)
    ny = Array(Int32, 1)
    ty = Array(Float64, nyest)
    c = Array(Float64, (nxest-kx-1) * (nyest-ky-1))
    fp = Array(Float64, 1)
    ier = Array(Int32, 1)
    
    # Work arrays.
    # Note that in lwrk, x and y are swapped with respect to the Fortran
    # documentation. See "NOTE REGARDING ARGUMENT ORDER" above. 
    lwrk = (4 + nyest * (mx+2*ky+5) + nxest * (2*kx+5) +
            my*(ky+1) + mx*(kx+1) + max(mx, nyest))
    wrk = Array(Float64, lwrk)
    kwrk = 3 + mx + my + nxest + nyest
    iwrk = Array(Int32, kwrk)
            
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
          &0f0, &my, y, &mx, x, z, &yb, &ye, &xb, &xe, &ky, &kx, &@compat(Float64(s)),
          &nyest, &nxest, ny, ty, nx, tx, c, fp,
          wrk, &lwrk, iwrk, &kwrk, ier)
    
    if !(ier[1] == 0 || ier[1] == -1 || ier[1] == -2)
        error(_fit2d_messages[ier[1]])
    end
            
    # Resize output arrays to the size actually used.
    resize!(tx, nx[1])
    resize!(ty, ny[1])
    resize!(c, (nx[1] - kx - 1) * (ny[1] - ky - 1))
    
    return Spline2D(tx, ty, c, kx, ky, fp[1])
end


# Evaluate spline at individual points
function evaluate(spline::Spline2D, x::Vector{Float64},
                  y::Vector{Float64})
    m = length(x)
    assert(length(y) == m)
    
    ier = Array(Int32, 1)
    lwrk = spline.kx + spline.ky + 2
    wrk = Array(Float64, lwrk)
    z = Array(Float64, m)
    
    ccall((:bispeu_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # ty, ny
           Ptr{Float64}, Ptr{Int32},  # tx, nx
           Ptr{Float64},  # c
           Ptr{Int32}, Ptr{Int32},  # ky, kx
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # y, x, z
           Ptr{Int32},  # m
           Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),  # wrk, lwrk, ier
          spline.ty, &length(spline.ty),
          spline.tx, &length(spline.tx),
          spline.c, &spline.ky, &spline.kx, y, x, z, &m,
          wrk, &lwrk, ier)
    
    ier[1] == 0 || error(_eval2d_message)
    
    return z
end

function evaluate(spline::Spline2D, x::Float64, y::Float64)
    ier = Array(Int32, 1)
    lwrk = spline.kx + spline.ky + 2
    wrk = Array(Float64, lwrk)
    z = Array(Float64, 1)
    ccall((:bispeu_, ddierckx), Void,
          (Ptr{Float64}, Ptr{Int32},  # ty, ny
           Ptr{Float64}, Ptr{Int32},  # tx, nx
           Ptr{Float64},  # c
           Ptr{Int32}, Ptr{Int32},  # ky, kx
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # y, x, z
           Ptr{Int32},  # m
           Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),  # wrk, lwrk, ier
          spline.ty, &length(spline.ty),
          spline.tx, &length(spline.tx),
          spline.c, &spline.ky, &spline.kx, &y, &x, z, &1,
          wrk, &lwrk, ier)
    
    ier[1] == 0 || error(_eval2d_message)
    return z[1]
end

# Evaluate spline on the grid spanned by the input arrays.
function evalgrid(spline::Spline2D, x::Vector{Float64},
                  y::Vector{Float64})
    mx = length(x)
    my = length(y)

    lwrk = mx*(spline.kx + 1) + my*(spline.ky + 1)
    wrk = Array(Float64, lwrk)
    kwrk = mx + my
    iwrk = Array(Int32, kwrk)
    ier = Array(Int32, 1)
    z = Array(Float64, mx, my)
    
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
          spline.ty, &length(spline.ty),
          spline.tx, &length(spline.tx),
          spline.c, &spline.ky, &spline.kx, y, &my, x, &mx, z,
          wrk, &lwrk, iwrk, &kwrk, ier)

    ier[1] == 0 || error(_eval2d_message)
    
    return z
end


end # module
