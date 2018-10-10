#!/usr/bin/env julia
using Dierckx
using Test
using Random: seed!

# Answers 'ans' are from scipy.interpolate,
# generated with genanswers.py script.

# -----------------------------------------------------------------------------
# Spline1D

x = [1., 2., 3.]
y = [0., 2., 4.]
spl = Spline1D(x, y; k=1, s=length(x))

yi = evaluate(spl, [1.0, 1.5, 2.0])
@test yi ≈ [0.0, 1.0, 2.0]
@test evaluate(spl, 1.5) ≈ 1.0
@test get_knots(spl) ≈ [1., 3.]
@test get_coeffs(spl) ≈ [0., 4.]
@test isapprox(get_residual(spl), 0.0, atol=1.e-30)

@test spl([1.0, 1.5, 2.0]) ≈ [0.0, 1.0, 2.0]
@test spl(1.5) ≈ 1.0

# test that a copy is returned by get_knots()
knots = get_knots(spl)
knots[1] = 1000.
@test get_knots(spl) ≈ [1., 3.]

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
spl = Spline1D(x, y; w=w, s=Float64(length(x)))
desired = [0.35100374, 0.51715855, 0.87789547, 0.98719344]
actual = evaluate(spl, [0.1, 0.5, 0.9, 0.99])
@test isapprox(actual, desired, atol=5e-4)

# test periodic
x = [1., 2., 3., 4., 5.]
y = [4., 1., 4., 1., 4.]
spl = Spline1D(x, y, periodic=true)

@test derivative(spl, 1) ≈ derivative(spl, 5)
@test derivative(spl, 1, nu=2) ≈ derivative(spl, 5, nu=2)

# tests for out-of-range
x = [0.0:4.0;]
y = x.^3

xp = range(-8.0, stop=13.0, length=100)
xp_zeros = Float64[(0. <= xi <= 4.) ? xi : 0.0 for xi in xp]
xp_clip = Float64[(0. <= xi <= 4.) ? xi : (xi < 0.0) ? 0.0 : 4. for xi in xp]

spl = Spline1D(x, y)
t = get_knots(spl)[2: end-1]  # knots, excluding those at endpoints
spl2 = Spline1D(x, y, t)

@test evaluate(spl, xp) ≈ xp_clip.^3
@test evaluate(spl2, xp) ≈ xp_clip.^3

# test other bc's
spl = Spline1D(x, y; bc="extrapolate")
@test evaluate(spl, xp) ≈ xp.^3
spl = Spline1D(x, y; bc="zero")
@test evaluate(spl, xp) ≈ xp_zeros.^3
spl = Spline1D(x, y; bc="error")
@test_throws ErrorException evaluate(spl, xp)

# test unknown bc
@test_throws ErrorException Spline1D(x, y; bc="unknown")

# test derivative
x = range(0, stop=1, length = 70)
y = x.^3
spl = Spline1D(x, y)
xt = [0.3, 0.4, 0.5]
@test derivative(spl, xt) ≈ 3xt.^2

# test integral
x = range(0, stop=10, length = 70)
y = x.^2
spl = Spline1D(x, y)
@test integrate(spl, 1.0, 5.0) ≈ 5.0^3/3 - 1/3

# test roots
x = range(0, stop=10, length = 70)
y = (x .- 4).^2 .- 1
spl = Spline1D(x, y)
@test roots(spl) ≈ [3, 5]

# test that show works.
io = IOBuffer()
show(io, spl)
seek(io, 0)
s = read(io, String)
@test s[1:9] == "Spline1D("

# -----------------------------------------------------------------------------
# ParametricSpline

u = [1., 2., 3.]
x = [1. 2. 3.; 0. 2. 4.]
spl = ParametricSpline(u, x, k=1, s=size(x, 2))

xi = evaluate(spl, [1.0, 1.5, 2.0])
@test xi ≈ [1.0 1.5 2.0; 0.0 1.0 2.0]
@test evaluate(spl, 1.5) ≈ [1.5, 1.0]
@test get_knots(spl) ≈ [1., 3.]
@test get_coeffs(spl) ≈ [1.0 3.0; 0.0 4.0]
@test isapprox(get_residual(spl), 0.0, atol=1.e-30)

@test spl([1.0, 1.5, 2.0]) ≈ [1.0 1.5 2.0; 0.0 1.0 2.0]
@test spl(1.5) ≈ [1.5, 1.0]

# test that a copy is returned by get_knots()
knots = get_knots(spl)
knots[1] = 1000.
@test get_knots(spl) ≈ [1., 3.]

# test periodic
x = [23. 24. 25. 25. 24. 23.;
     13. 12. 12. 13. 13. 13.]
spl = ParametricSpline(x, periodic=true)
@test evaluate(spl, 0) ≈ evaluate(spl, 1)
@test derivative(spl, 0) ≈ derivative(spl, 1)
@test derivative(spl, 0, nu=2) ≈ derivative(spl, 1, nu=2)

# tests for out-of-range
u = 0.0:4.0
x = [u'.^2; u'.^3]

up = range(-8.0, stop=13.0, length = 100)
up_zeros = Float64[(0. <= ui <= 4.) ? ui : 0.0 for ui in up]
up_clip = Float64[(0. <= ui <= 4.) ? ui : (ui < 0.0) ? 0.0 : 4. for ui in up]

spl = ParametricSpline(u, x)
t = get_knots(spl)[2: end-1]  # knots, excluding those at endpoints
spl2 = ParametricSpline(u, x, t)

@test evaluate(spl, up) ≈ [up_clip'.^2; up_clip'.^3]
@test evaluate(spl2, up) ≈ [up_clip'.^2; up_clip'.^3]

# test other bc's
spl = ParametricSpline(u, x; bc="extrapolate")
@test evaluate(spl, up) ≈ [up'.^2; up'.^3]
spl = ParametricSpline(u, x; bc="zero")
@test evaluate(spl, up) ≈ [up_zeros'.^2; up_zeros'.^3]
spl = ParametricSpline(u, x; bc="error")
@test_throws ErrorException evaluate(spl, up)

# test unknown bc
@test_throws ErrorException ParametricSpline(u, x; bc="unknown")

# test derivative
u = range(0, stop=1, length = 70)
x = [u'.^2; u'.^3]
spl = ParametricSpline(u, x)
ut = [0.3, 0.4, 0.5]
@test derivative(spl, 0.3) ≈ [2*0.3, 3*0.3^2]
@test derivative(spl, ut) ≈ [2*ut'; 3*ut'.^2]
@test derivative(spl, 0.3, nu=2) ≈ [2.0, 6*0.3]
@test derivative(spl, ut, nu=2) ≈ [2*ones(3)'; 6*ut']

# test integral
u = range(0, stop=10, length = 70)
x = [u'.^2; u'.^3]
spl = ParametricSpline(u, x)
@test integrate(spl, 1.0, 5.0) ≈ [5.0^3/3 - 1/3, 5.0^4/4 - 1/4]

# test that show works.
io = IOBuffer()
show(io, spl)
seek(io, 0)
s = read(io, String)
@test s[1:17] == "ParametricSpline("

# -----------------------------------------------------------------------------
# Spline2D

# test linear
x = [1., 1., 1., 2., 2., 2., 3., 3., 3.]
y = [1., 2., 3., 1., 2., 3., 1., 2., 3.]
z = [0., 0., 0., 2., 2., 2., 4., 4., 4.]
spl = Spline2D(x, y, z; kx=1, ky=1, s=length(x))
tx, ty = get_knots(spl)
@test tx ≈ [1., 3.]
@test ty ≈ [1., 3.]
@test isapprox(get_residual(spl), 0.0, atol=1e-16)
@test evaluate(spl, 2.0, 1.5) ≈ 2.0
@test evalgrid(spl, [1.,1.5,2.], [1.,1.5]) ≈ [0. 0.; 1. 1.; 2. 2.]

# test 1-d grid arrays
@test evalgrid(spl, [2.0], [1.5])[1, 1] ≈ 2.0

# In this setting, lwrk2 is too small in the default run.
x = range(-2, stop=2, length = 80)
y = range(-2, stop=2, length = 80)
z = x .+ y
spl = Spline2D(x, y, z; s=length(x))
@test evaluate(spl, 1.0, 1.0) ≈ 2.0

# In this setting lwrk2 is too small multiple times!
# Eventually an error about s being too small is thrown.
seed!(0)
x = rand(100)
y = rand(100)
z = sin.(x) .* sin.(y)
@test_throws ErrorException Spline2D(x, y, z; kx=1, ky=1, s=0.0)

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
@test zi ≈ ans

zi = spl(xi, yi)
@test zi ≈ ans

# grid output
xi = [1., 1.5, 2.3, 4.5]
yi = [1., 2.3, 5.3]
ans = [2.94429906542  1.16946130841  1.99831775701;
       2.80393858478  1.25537598131  1.99873831776;
       1.67143209613  1.94853338542  2.00063588785;
       1.89392523364  1.8126946729  2.01042056075]
zi = evalgrid(spl, xi, yi)
@test zi ≈ ans


# Test 2-d integration
test2d_1(x, y) = 1 - x^2 -y^2
test2d_2(x, y) = cos(x) + sin(y)
test2d_3(x, y) = x*exp(x-y)
for (f, domain, exact) in [(test2d_1, (0.0, 1.0, 0.0, 1.0), 1.0/3.0),
                           (test2d_2, (0.0, pi, 0.0, pi), 2.0*pi),
                           (test2d_3, (0.0, 1.0, 0.0, 1.0), (ℯ-1.0)/ℯ)]
    (x0, x1, y0, y1) = domain

    # define grids for x and y dimensions:
    npoints = 50
    xgrid = range(x0, stop=x1, length = npoints)
    ygrid = range(y0, stop=y1, length = npoints)

    fxygrid = zeros(npoints, npoints)
    for (j, y) in enumerate(ygrid)
        local j, y
        for (i, x) in enumerate(xgrid)
            local i, x
            fxygrid[i,j] = f(x, y)
        end
    end

    spl1 = Spline2D(xgrid, ygrid, fxygrid)
    @test isapprox(integrate(spl1, x0, x1, y0, y1), exact, atol=1e-6)
end

println("All tests passed.")
