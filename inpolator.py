__author__ = 'p0054421'

from scipy.interpolate import RegularGridInterpolator
import numpy as np
from scipy.interpolate import griddata


def interp(nx, ny, nz, velp, dir, t, rr):
    x = np.linspace(0, nx - 1, nx)
    y = np.linspace(0, ny - 1, ny)
    z = np.linspace(0, nz - 1, nz)

    print x.shape, 'xshape'

    # ) * sin(Y*2*pi*3) *
    # a = np.exp(-0.5 * (((X-0.3)/0.2)**2 + ((Y-0)/0.2)**2 * 0 + ((Z-0.6)/0.3)**2)) * (sin(Y*2*pi*3) + 1)
    a = velp[:, :, :, dir, t]
    # print a.shape, a.min(), a.max()
    xv = np.linspace(0, nx - 1, rr * nx)
    yv = np.linspace(0, ny - 1, rr * ny)
    zv = np.linspace(0, nz - 1, rr * nz)
    toto = np.swapaxes(np.array([xv, yv, zv]), 1, 0)
    # totoi= np.swapaxes(np.array([x, y, z]),1,0)
    # toto=toto.ravel()
    k = 0
    toto = np.empty((rr * rr * nx * ny, 2))
    for i in xv:
        for j in yv:
            toto[k, :] = i, j
            k += 1

    print toto.shape, xv.shape
    # xv = x
    # yv = y
    # zv = z
    print 'evaluating'
    # np.array([totoi[0],totoi[1],totoi[2]]).shape

    # my_interpolating_function = RegularGridInterpolator((x, y, z), a, method='linear', bounds_error=False, fill_value=0.0)
    # my_interpolating_function1 = griddata(totoi, a[0], toto, method='linear')
    # print a.shape, my_interpolating_function1.shape, toto.shape
    print 'end'

    return my_interpolating_function1
