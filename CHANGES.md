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
