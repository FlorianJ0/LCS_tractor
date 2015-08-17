__author__ = 'p0054421'

from scipy.interpolate import interpn
import numpy as np
from numpy import sin, pi

def interp(nx,ny,nz,velp,dir,t):
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    z = np.linspace(0, nz-1, nz)

    print x.shape, 'xshape'


    Z, Y = np.meshgrid(y, z)
    X = np.array([[x]]).transpose()
    X = Y * 0 + X
    Y = X * 0 + Y
    Z = X * 0 + Z
    print X.shape, Y.shape, Z.shape
    # ) * sin(Y*2*pi*3) *
    a = np.exp(-0.5 * (((X-0.3)/0.2)**2 + ((Y-0)/0.2)**2 * 0 + ((Z-0.6)/0.3)**2)) * (sin(Y*2*pi*3) + 1)
    a = velp[:,:,:,dir,t]
    print a.shape, a.min(), a.max()
    # xv = np.linspace(0, 1, 33)
    # yv = np.linspace(0, 1, 33)
    # zv = np.linspace(0, 1, 33)
    xv = x
    yv = y
    zv = z

    print 'evaluating'
    values = []
    for xi in xv:
        yvalues = []
        for yi in yv:
            zvalues = []
            for zi in zv:
                v = interpn(x, y, z, a, xi, zi, yi)
                #z = np.exp(-0.5 * (((xi-0.3)/0.2)**2 + ((yi-0.6)/0.3)**2))
                zvalues.append(v)
            yvalues.append(zvalues)
        values.append(yvalues)
    values = np.array(values)
    print 'plotting'
