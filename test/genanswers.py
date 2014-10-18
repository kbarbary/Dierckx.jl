#!/usr/bin/env python
"""generate 'correct' answers with scipy"""
import numpy as np
from scipy.interpolate import RectBivariateSpline

def print_array_jl_syntax(arr):
    """Print an array to stdout in julia syntax, copy-pastable into a
    julia script"""

    if arr.ndim == 1:
        ending = ","
        nx = arr.shape[0]
        print "ans = [" + str(arr[0]) + ending
        for i in range(1, nx):
            if i == nx - 1:
                ending = "]"
            print "       " + str(arr[i]) + ending

    elif arr.ndim == 2:
        ending = ";"
        nx, ny = arr.shape
        print "ans = [" + "  ".join([str(arr[0, j]) for j in range(ny)]) + ending
        for i in range(1, nx):
            if i == nx - 1:
                ending = "]"
            print "       " + "  ".join([str(arr[i, j]) for j in range(ny)]) + ending



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

# element-wise output
xi = np.array([1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.])
yi = np.array([1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.])
zi = s(xi, yi, grid=False)
print_array_jl_syntax(zi)

# Grid output
print "\n"
xi = np.array([1., 1.5, 2.3, 4.5])
yi = np.array([1., 2.3, 5.3])
zi = s(xi, yi, grid=True)
print_array_jl_syntax(zi)
