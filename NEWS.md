v0.4.1 (2018-10-10)
===================

- Fix bug when displaying Spline1D and ParametricSpline on Julia 1.0.
- Use gcc multiversioning in BinaryBuilder to get a
  GCC-version-specific libddierckx binary. This addresses problems
  where libddireckx was built against libgfortran4 (which is
  distributed with the official Julia 1.0 binaries, built with GCC 7),
  but users building Julia from source with GCC 8 had libgfortran5.

v0.4.0 (2018-08-12)
===================

- Drop support for Julia 0.6 and before (support Julia 0.7+ only)
- Use pre-compiled binaries on all supported platforms.

v0.3.0 (2017-09-18)
===================

- Drop support for Julia 0.4
- Fix most Julia 0.7 deprecation warnings

v0.2.1 (2017-06-01)
===================

- Add `ParametricSpline` type.

v0.2.0 (2017-01-26)
===================

- Fix Julia 0.6 deprecation warnings
- Drop support for Julia 0.3

v0.1.9 (2016-08-10)
===================

- Fix Julia v0.5 deprecation warnings
- Add flexibility in Fortran compiler (UNIX only) in build script

v0.1.8 (2016-02-29)
===================

- Enable precompilation on v0.4
- Use BinDeps download command on Windows in build script.

v0.1.7 (2015-09-28)
===================

- Fix v0.4 deprecation warnings
- Add `call` methods as synonyms for `evaluate()` on v0.4+,
  allowing syntax such as `spline(x)`.

v0.1.6 (2015-07-27)
===================

- New download location for windows binaries.

v0.1.5 (2015-06-08)
===================

- Adds `integrate` and `roots` methods for `Spline1D`.
- Widen allowed argument types to `AbstractArray` in spline construction
  and evaluation. Internally, arguments are converted to Float64 arrays (if
  needed).

v0.1.4 (2015-03-23)
===================

- Fix v0.4 deprecation warnings
- Add `__init__` method to enable precompilation.

v0.1.3 (2014-12-01)
===================

- Fix bug in Spline2D(::Vector, ::Vector, ::Vector) where an error
  code from the Fortran library indicating that the second work array is too
  small was not being handled.
- add show() method for Spline1D
- add derivative() method for Spline1D

v0.1.2 (2014-11-07)
===================

- Add support for Windows by downloading a compiled dll.

v0.1.1 (2014-10-27)
===================

- Fix Makefile bug in OS X build.

v0.1.0 (2014-10-24)
===================

Initial release.
