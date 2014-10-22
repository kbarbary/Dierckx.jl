#!/usr/bin/env julia
using Dierckx
using Base.Test

# Answers 'ans' are from scipy.interpolate,
# generated with genanswers.py script.

x = [0.5, 2., 3., 4., 5.5, 8.]
y = [0.5, 2., 3., 4.]

# shape is (nx, ny)
z = [1. 2. 1. 2.;
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
