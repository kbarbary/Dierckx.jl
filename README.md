Dierckx.jl
==========

*Julia library for 1-d and 2-d splines*

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
julia> Pkg.add("git://github.com/kbarbary/Dierckx.jl.git")
```

Example Usage
-------------

```julia
using Dierckx

x = [0.5, 2., 3., 4., 5.5, 8.]
y = [0.5, 2., 3., 4.]

# size of z is (length(x), length(y))
z = [1. 2. 1. 2.;
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
Spline1D(x, y; w=ones(length(x)), k=3, s=0.0)
Spline1D(x, y, xknots; w=ones(length(x)), k=3)
```

Create a spline from vectors `x` and `y` of degree `k` (1 = linear, 2
= quadratic, 3 = cubic, etc). In the first form, the position and
number of knots is chosen automatically. The smoothness of the spline
is then achieved by minimalizing the discontinuity jumps of the `k`th
derivative of the spline at the knots. The amount of smoothness is
determined by the condition that `sum((w[i]*(y[i]-spline(x[i])))**2)
<= s`, with `s` a given non-negative constant, called the smoothing
factor. By means of this parameter, the user can control the tradeoff
between closeness of fit and smoothness of fit of the approximation.
if `s` is too large, the spline will be too smooth and signal will be
lost ; if `s` is too small the spline will pick up too much noise. in
the extreme cases the program will return an interpolating spline if
`s=0.0` and the weighted least-squares polynomial of degree `k` if `s`
is very large.

In the second form, the knots are supplied by the user. There is
no smoothing parameter in this form.

```julia
evaluate(spl, x)
```
Evalute the 1-d spline `spl` at points given in `x`, which can be a
Vector or scalar. If a Vector, the input arrays must be monotonically
increasing.

### 2-d Splines

```julia
Spline2D(x, y, z; kx=3, ky=3, s=0.0)
```

Fit a 2-d spline to data points in `x`, `y`, `z`. `x` and `y` must be
vectors.

If `z` is also a vector, the inputs are assumed to represent
unstructured data, with `z[i]` being the function value at point
`(x[i], y[i])`. In this case, the lengths of all inputs must match.

If `z` is a 2-d array, the data are assumed to be gridded: `z[i, j]`
is the function value at `(x[i], y[j])`. In this case, it is required
that `size(z) == (length(x), length(y))`.

```julia
evaluate(spl, x, y)
```

Evalute the 2-d spline `spl` at points `(x[i], y[i])`. Points outside
the domain of the spline are set to the values at the boundary.

```julia
evalgrid(spl, x, y)
```

Evaluate the 2-d spline `spl` at the grid points spanned by the
coordinate arrays `x` and `y`. The input arrays must be monotonically
increasing.

Translation from scipy.interpolate
----------------------------------

The `*Spline` classes in scipy.interpolate are also thin wrappers
for the Dierckx Fortran library. The performance of Dierckx.jl should
be similar or better than the `scipy.interpolate` classes. (Better for
small arrays where Python overhead is more significant.) The
equivalent of a specific classes in `scipy.interpolate`:

| scipy.interpolate class      | Dierckx.jl constructor method            |
| ---------------------------- | ---------------------------------------- |
| UnivariateSpline             | Spline1D(x, y)                           |
| InterpolatedUnivariateSpline | Spline1D(x, y; s=0.0) [s=0.0 is default] |
| LSQUnivariateSpline          | Spline1D(x, y, xknots)                   |
| SmoothBivariateSpline        | Spline2D()                               |
| LSQBivariateSpline           |                                          |
| RectBivariateSpline          | Spline2D() with 2-d z array              |
| SmoothSphereBivariateSpline  |                                          |
| LSQSphereBivariateSpline     |                                          |
| RectSphereBivariateSpline    |                                          |



License
-------

Dierckx.jl is distributed under a 3-clause BSD license. See LICENSE.md
for details. The real*8 version of the Dierckx Fortran library as well as
some test cases and error messages are copied from the scipy package,
which is distributed under this license.
