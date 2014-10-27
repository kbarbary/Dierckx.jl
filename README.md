Dierckx.jl
==========

*Julia library for 1-d and 2-d splines*

[![Build Status](https://travis-ci.org/kbarbary/Dierckx.jl.svg?branch=master)](https://travis-ci.org/kbarbary/Dierckx.jl)
[![Coverage Status](https://coveralls.io/repos/kbarbary/Dierckx.jl/badge.png?branch=master)](https://coveralls.io/r/kbarbary/Dierckx.jl?branch=master)


This is a Julia wrapper for the
[dierckx](http://www.netlib.org/dierckx/index.html) Fortran library,
the same library underlying the spline classes in scipy.interpolate.
Some of the functionality here overlaps with [Grid.jl](
https://github.com/timholy/Grid.jl), a pure-Julia interpolation
package. This package is intended to complement Grid.jl and to serve
as a benchmark in cases of overlapping functionality.

Install
-------

```julia
julia> Pkg.add("Dierckx")
```

The Fortran library source code is distributed with the package, so
you need a Fortran compiler. On Ubuntu, `sudo apt-get install gfortran`
will do it.

Windows is not yet supported.

Example Usage
-------------

```julia
using Dierckx

# -------------------------------------------------------
# 1-d splines

x = [0., 1., 2., 3., 4.]
y = [0., 1., 8., 27., 64.]
spl = Spline1D(x, y)
evaluate(spl, [1.5, 2.5])  # result = [3.375, 15.625]
evaluate(spl, 1.5)  # result = 3.375

# -------------------------------------------------------
# 2-d splines

x = [0.5, 2., 3., 4., 5.5, 8.]
y = [0.5, 2., 3., 4.]
z = [1. 2. 1. 2.;  # size is (length(x), length(y))
     1. 2. 1. 2.;
     1. 2. 3. 2.;
     1. 2. 2. 2.;
     1. 2. 1. 2.;
     1. 2. 3. 1.]

# Fit a 2-d spline to inputs
spline = Spline2D(x, y, z)

# evaluate at element-wise points
xi = [1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.]
yi = [1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.]
zi = evaluate(spline, xi, yi)  # 1-d array of length 7

# evaluate at grid spanned by input arrays
xi = [1., 1.5, 2.3, 4.5]
yi = [1., 2.3, 5.3]
zi = evalgrid(spline, xi, yi)  # 2-d array of size (4, 3)
```

Reference
---------

### 1-d Splines

```julia
Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
Spline1D(x, y, xknots; w=ones(length(x)), k=3, bc="nearest")
```

Create a spline of degree `k` (1 = linear, 2 = quadratic, 3 = cubic, up to 5)
from 1-d arrays `x` and `y`. The parameter `bc` specifies the behavior when
evaluating the spline outside the support domain, which is `(minimum(x),
maximum(x))`. The allowed values are `"nearest"`, `"zero"`, `"extrapolate"`,
`"error"`.

In the first form, the number and positions of knots are chosen
automatically. The smoothness of the spline is then achieved by
minimalizing the discontinuity jumps of the `k`th derivative of the
spline at the knots. The amount of smoothness is determined by the
condition that `sum((w[i]*(y[i]-spline(x[i])))**2) <= s`, with `s` a
given non-negative constant, called the smoothing factor. The number
of knots is increased until the condition is satisfied. By means of
this parameter, the user can control the tradeoff between closeness of
fit and smoothness of fit of the approximation.  if `s` is too large,
the spline will be too smooth and signal will be lost ; if `s` is too
small the spline will pick up too much noise. in the extreme cases the
program will return an interpolating spline if `s=0.0` and the
weighted least-squares polynomial of degree `k` if `s` is very large.

In the second form, the knots are supplied by the user. There is no
smoothing parameter in this form. The program simply minimizes the
discontinuity jumps of the `k`th derivative of the spline at the given
knots.

```julia
evaluate(spl, x)
```
Evalute the 1-d spline `spl` at points given in `x`, which can be a
1-d array or scalar. If a 1-d array, the values must be monotonically
increasing.

### 2-d Splines

```julia
Spline2D(x, y, z; w=ones(length(x)), kx=3, ky=3, s=0.0)
Spline2D(x, y, z; kx=3, ky=3, s=0.0)
```

Fit a 2-d spline to the input data. `x` and `y` must be 1-d arrays.

If `z` is also a 1-d array, the inputs are assumed to represent
unstructured data, with `z[i]` being the function value at point
`(x[i], y[i])`. In this case, the lengths of all inputs must match.

If `z` is a 2-d array, the data are assumed to be gridded: `z[i, j]`
is the function value at `(x[i], y[j])`. In this case, it is required
that `size(z) == (length(x), length(y))`.

```julia
evaluate(spl, x, y)
```

Evalute the 2-d spline `spl` at points `(x[i], y[i])`. Inputs can be
Vectors or scalars. Points outside the domain of the spline are set to
the values at the boundary.

```julia
evalgrid(spl, x, y)
```

Evaluate the 2-d spline `spl` at the grid points spanned by the
coordinate arrays `x` and `y`. The input arrays must be monotonically
increasing.

Translation from scipy.interpolate
----------------------------------

The `Spline` classes in scipy.interpolate are also thin wrappers
for the Dierckx Fortran library. The performance of Dierckx.jl should
be similar or better than the scipy.interpolate classes. (Better for
small arrays where Python overhead is more significant.) The
equivalent of a specific classes in scipy.interpolate:

| scipy.interpolate class      | Dierckx.jl constructor method              |
| ---------------------------- | ------------------------------------------ |
| UnivariateSpline             | `Spline1D(x, y; s=length(x))`              |
| InterpolatedUnivariateSpline | `Spline1D(x, y; s=0.0)`                    |
| LSQUnivariateSpline          | `Spline1D(x, y, xknots)`                   |
| SmoothBivariateSpline        | `Spline2D(x, y, z; s=length(x))`           |
| LSQBivariateSpline           |                                            |
| RectBivariateSpline          | `Spline2D(x, y, z; s=0.0)` (z = 2-d array) |
| SmoothSphereBivariateSpline  |                                            |
| LSQSphereBivariateSpline     |                                            |
| RectSphereBivariateSpline    |                                            |



License
-------

Dierckx.jl is distributed under a 3-clause BSD license. See LICENSE.md
for details. The real*8 version of the Dierckx Fortran library as well as
some test cases and error messages are copied from the scipy package,
which is distributed under this license.
