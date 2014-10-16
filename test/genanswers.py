#!/usr/bin/env python
"""generate 'correct' answers with scipy"""
import numpy as np
from scipy.interpolate import RectBivariateSpline

x = np.array([0.5, 2., 3., 4., 5.5, 8.])
y = np.array([0.5, 2., 3., 4.])

# z has shape (nx, ny) in *python*
z = np.array([[1., 2., 1., 2.],
              [1., 2., 1., 2.],
              [1., 2., 3., 2.],
              [1., 2., 2., 2.],
              [1., 2., 1., 2.],
              [1., 2., 3., 1.]])

s = RectBivariateSpline(x, y, z)

xi = np.array([1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.])
yi = np.array([1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.])

zi = s(xi, yi, grid=False)


print "ans = [" + str(zi[0]) + ","
ending = ","
for i in range(1, len(zi)):
    if i == len(zi) - 1:
        ending = "]"
    print "       " + str(zi[i]) + ending
