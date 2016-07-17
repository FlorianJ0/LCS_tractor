import numpy as np

from interpolation.splines import CubicSpline

a = np.array([0.0,0.0,0.0])         # lower boundaries
b = np.array([1.0,1.0,1.0])         # upper boundaries
orders = np.array([50,50,50])       # 10 points along each dimension
values = np.random.random(orders)   # values at each node of the grid
# print values.shape, v0.shape
S = np.random.random((10^6,3))    # coordinates at which to evaluate the splines

# # multilinear
# lin = LinearSpline(a,b,orders,values)
# #
# V = lin(S)
# print V
# # quit()

# cubic
spline = CubicSpline(a,b,orders,values) # filter the coefficients
V = spline(S)
print V# interpolates -> (100000,) array

