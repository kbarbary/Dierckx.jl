#!/usr/bin/env julia
using Dierckx
using Base.Test


x = [0.5, 2., 3., 4., 5.5, 8.]
y = [0.5, 2., 3., 4.]

# shape is (nx, ny)
z = [1. 2. 1. 2.;
     1. 2. 1. 2.;
     1. 2. 3. 2.;
     1. 2. 2. 2.;
     1. 2. 1. 2.;
     1. 2. 3. 1.]

s = GridSpline(x, y, z)

xi = [1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.]
yi = [1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.]

zi = evaluate(s, xi, yi)

# Answer from scipy.interpolate.RectBivariateSpline,
# generated with genanswers
ans = [2.94429906542,
       1.25537598131,
       2.00063588785,
       1.0,
       2.93952664,
       1.06482509358,
       3.0]

@test_approx_eq(zi, ans)
