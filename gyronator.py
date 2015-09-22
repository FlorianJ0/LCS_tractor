__author__ = 'p0054421'

import numpy as np

nx = 50
ny = 25
nz = 25
dx = 4e-2
dy = 4e-2
dz = 1e-2
ttot = 5
dt = 0.5
w = 2 * np.pi / 10


# w=0
def a(t):
    r = 0.5 * np.sin(w * t)
    return r


def b(t):
    s = 1 - 2 * 0.5 * np.sin(w * t)
    return s


def f(x, t):
    toto = a(t) * x ** 2 + b(t) * x
    return toto


def dfdx(x, t):
    toto = 2 * a(t) * x + b(t)
    return toto


def gyro():
    vel = np.empty((nx, ny, 2, int(ttot / dt)))  #
    # nz = 5
    vel3D = np.zeros((nx, ny, nz, 3, int(ttot / dt))) * (-1e-3)  #

    for t in range(0, int(ttot / dt)):
        tt = t * dt
        print tt
        for i in xrange(nx):
            for j in xrange(ny):
                x = i * dx
                y = j * dy

                vel[i, j, 0, t] = - np.pi * 0.1 * np.sin(np.pi * f(x, tt)) * np.cos(np.pi * y) * 1
                vel[i, j, 1, t] = np.pi * 0.1 * np.cos(np.pi * f(x, tt)) * np.sin(np.pi * y) * dfdx(x, tt) * 1
        vel3D[:, :, 0, 0:2, t] = vel[:, :, :, t]

    for k in range(1, nz):
        vel3D[:, :, k, :, :] = vel3D[:, :, 0, :, :] -1e-3

    # vel3D[:, :, :, 2, :] = -1e-3

    print vel.shape
    domain = np.array([0, 2, 0, 1, 0, 0.1])
    # nz = 0
    dim_initial = 3
    print vel3D.shape, nx, ny, nz, domain, dim_initial
    return vel3D, nx, ny, nz, dim_initial, ttot, dt, domain
