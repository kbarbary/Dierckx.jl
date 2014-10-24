#!/usr/bin/env julia
using Dierckx
using Base.Test

# Answers 'ans' are from scipy.interpolate,
# generated with genanswers.py script.

# -----------------------------------------------------------------------------
# Spline1D

x = [1., 2., 3.]
y = [0., 2., 4.]
spl = Spline1D(x, y; k=1, s=length(x))

yi = evaluate(spl, [1.0, 1.5, 2.0])
@test_approx_eq(yi, [0.0, 1.0, 2.0])
@test_approx_eq(evaluate(spl, 1.5), 1.0)
@test_approx_eq(get_knots(spl), [1., 3.])
@test_approx_eq(get_coeffs(spl), [0., 4.])
@test_approx_eq(get_residual(spl), 0.0)

# test that a copy is returned by get_knots()
knots = get_knots(spl)
knots[1] = 1000.
@test_approx_eq(get_knots(spl), [1., 3.])

# test ported from scipy.interpolate testing this bug:
# http://mail.scipy.org/pipermail/scipy-dev/2008-March/008507.html
x = [-1., -0.65016502, -0.58856235, -0.26903553, -0.17370892,
     -0.10011001, 0., 0.10011001, 0.17370892, 0.26903553, 0.58856235,
     0.65016502, 1.]
y = [1.,0.62928599, 0.5797223, 0.39965815, 0.36322694, 0.3508061,
     0.35214793, 0.3508061, 0.36322694, 0.39965815, 0.5797223,
     0.62928599, 1.]
w = [1.00000000e+12, 6.88875973e+02, 4.89314737e+02, 4.26864807e+02,
     6.07746770e+02, 4.51341444e+02, 3.17480210e+02, 4.51341444e+02,
     6.07746770e+02, 4.26864807e+02, 4.89314737e+02, 6.88875973e+02,
     1.00000000e+12]
spl = Spline1D(x, y; w=w, s=float64(length(x)))
desired = [0.35100374, 0.51715855, 0.87789547, 0.98719344]
actual = evaluate(spl, [0.1, 0.5, 0.9, 0.99])
@test_approx_eq_eps(actual, desired, 5e-4)

# tests for out-of-range
x = [0.0:4.0]
y = x.^3

xp = linspace(-8.0, 13.0, 100)
xp_zeros = Float64[(0. <= xi <= 4.) ? xi : 0.0 for xi in xp]
xp_clip = Float64[(0. <= xi <= 4.) ? xi : (xi<0.0)? 0.0 : 4. for xi in xp]

spl = Spline1D(x, y)
t = get_knots(spl)[2: end-1]  # knots, excluding those at endpoints
spl2 = Spline1D(x, y, t)

@test_approx_eq(evaluate(spl, xp), xp_clip.^3)
@test_approx_eq(evaluate(spl2, xp), xp_clip.^3)

# test other bc's
spl = Spline1D(x, y; bc="extrapolate")
@test_approx_eq(evaluate(spl, xp), xp.^3)
spl = Spline1D(x, y; bc="zero")
@test_approx_eq(evaluate(spl, xp), xp_zeros.^3)
spl = Spline1D(x, y; bc="error")
@test_throws ErrorException evaluate(spl, xp)

# test unknown bc
@test_throws ErrorException Spline1D(x, y; bc="unknown")

# -----------------------------------------------------------------------------
# Spline2D

# test linear
x = [1., 1., 1., 2., 2., 2., 3., 3., 3.]
y = [1., 2., 3., 1., 2., 3., 1., 2., 3.]
z = [0., 0., 0., 2., 2., 2., 4., 4., 4.]
spl = Spline2D(x, y, z; kx=1, ky=1, s=length(x))
tx, ty = get_knots(spl)
@test_approx_eq(tx, [1., 3.])
@test_approx_eq(ty, [1., 3.])
@test_approx_eq_eps(get_residual(spl), 0.0, 1e-16)
@test_approx_eq(evaluate(spl, 2.0, 1.5), 2.0)
@test_approx_eq(evalgrid(spl, [1.,1.5,2.], [1.,1.5]), [0. 0.; 1. 1.; 2. 2.])

# test 1-d grid arrays
@test_approx_eq(evalgrid(spl, [2.0], [1.5])[1, 1], 2.0)

# In this setting, lwrk2 is too small in the default run.
x = linspace(-2, 2, 80)
y = linspace(-2, 2, 80)
z = x .+ y
spl = Spline2D(x, y, z; s=length(x))
@test_approx_eq(evaluate(spl, 1.0, 1.0), 2.0)

# test grid input creation
x = [0.5, 2., 3., 4., 5.5, 8.]
y = [0.5, 2., 3., 4.]
z = [1. 2. 1. 2.;  # shape is (nx, ny)
     1. 2. 1. 2.;
     1. 2. 3. 2.;
     1. 2. 2. 2.;
     1. 2. 1. 2.;
     1. 2. 3. 1.]
spl = Spline2D(x, y, z)

# element-wise output
xi = [1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.]
yi = [1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.]
ans = [2.94429906542,
       1.25537598131,
       2.00063588785,
       1.0,
       2.93952664,
       1.06482509358,
       3.0]
zi = evaluate(spl, xi, yi)
@test_approx_eq(zi, ans)

# grid output
xi = [1., 1.5, 2.3, 4.5]
yi = [1., 2.3, 5.3]
ans = [2.94429906542  1.16946130841  1.99831775701;
       2.80393858478  1.25537598131  1.99873831776;
       1.67143209613  1.94853338542  2.00063588785;
       1.89392523364  1.8126946729  2.01042056075]
zi = evalgrid(spl, xi, yi)
@test_approx_eq(zi, ans)

println("All tests passed.")
