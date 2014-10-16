Dierckx.jl
==========

A Julia wrapper for the
[Dierckx](http://www.netlib.org/dierckx/index.html) Fortran library
for spline fitting.

The functionality in the Dierckx library overlaps with [Grid.jl](
https://github.com/timholy/Grid.jl), a pure-Julia interpolation
package. This package is intended to be used as a benchmark (in both
performance and accuracy), and to complement Grid.jl where Dierckx has
some unique functionality.

Note: `scipy.interpolate` is also a thin wrapper for the same Fortran
Dierckx library. The functionality and performance here should be
similar. Equivalent classes/types:

| scipy.interpolate   | Dierckx.jl |
| ------------------- | ---------- |
| RectBivariateSpline | GridSpline |

_Note that Dierckx.jl is in development and type names may change._

Usage
-----

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

spline = GridSpline(x, y, z)

xi = [1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.]
yi = [1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.]

zi = evaluate(spline, xi, yi)  # 1-d array of length 7
```

License
-------

Dierckx.jl is distributed under a 3-clause BSD license. See LICENSE.md
for details.

