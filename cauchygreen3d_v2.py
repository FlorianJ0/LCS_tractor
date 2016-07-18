__author__ = 'p0054421'
__date__ = '30 Juillet 2015'

import ConfigParser

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
from scipy import *

from interpolation.splines import LinearSpline, CubicSpline

# noinspection PyPep8Naming
Config = ConfigParser.ConfigParser()
Config.read('parameters.ini')
# cubic interp ? if 0 then linear ADVECTION
cubicADV = 1
# cubic interp ? if 0 then linear CGtensor
cubicCG = 0


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def rpog():
    print("# "),


def cgstki3(velp, zplan, tt, dt, nx, ny, nz, dim, domain, simtstep):
    # x = np.arange(tt)
    print 'domain :'
    print domain
    ttt = int(tt / dt)
    ttt = np.arange(ttt)
    tmin = ttt[0] * dt
    tmax = ttt[-1] * dt
    print 'tmin, tmax', tmin, tmax
    # print ttt[0], ttt[-1]
    nt = len(ttt)
    print 'n time step', nt

    # small (so is your dick) vector d1(d1 0 0) d2(0 d2 0) d3(0 0 d3)
    minx, maxx, miny, maxy, minz, maxz = domain

    dx = (maxx - minx) / nx
    dy = (maxy - miny) / ny
    dz = (maxz - minz) / nz
    print 'dx=', dx
    if abs(dx - dy) > 1e-5 or abs(dx - dz) > 1e-5:
        print 'bad spacing', dx, dy, dz
        quit()
    a = np.array([minx, miny, minz, tmin])  # lower boundaries
    b = np.array([maxx, maxy, maxz, tmax])  # upper boundaries

    orders = np.array([nx, ny, nz, nt])  # 10 points along each dimension

    velpu = velp[:, :, :, 0, :]
    velpv = velp[:, :, :, 1, :]
    velpw = velp[:, :, :, 2, :]
    if cubicADV:
        linx = CubicSpline(a, b, orders, velpu)
        liny = CubicSpline(a, b, orders, velpv)
        linz = CubicSpline(a, b, orders, velpw)
    else:
        linx = LinearSpline(a, b, orders, velpu)
        liny = LinearSpline(a, b, orders, velpv)
        linz = LinearSpline(a, b, orders, velpw)
    # intitial pos:
    grid_x, grid_y, grid_z, grid_t = np.mgrid[minx:maxx - dx:nx * 1j, miny:maxy - dx:ny * 1j,
                                     minz:maxz - dx:nz * 1j, tmin:tmax:nt * 1j]

    positions = np.vstack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel(), grid_t.ravel()]).T

    # interpUZ = linz(positions)
    # interpUX = linx(positions)
    # interpUY = liny(positions)

    pos_x, pos_y, pos_z = np.mgrid[minx:maxx - dx:nx * 1j, miny:maxy - dx:ny * 1j,
                          minz:maxz - dx:nz * 1j]
    pos_x = pos_x.ravel()
    pos_y = pos_y.ravel()
    pos_z = pos_z.ravel()

    # cheap advection

    t0 = tmin
    z0 = positions
    t1 = tmax
    dt = 0.005  # interp step
    listt = np.linspace(t0, t1, int((tmax - tmin) / dt))
    print 'number of interp steps', len(listt)
    for t in listt:
        aa = np.empty([pos_x.shape[0], 4])
        aa[:, 0] = pos_x
        aa[:, 1] = pos_y
        aa[:, 2] = pos_z
        aa[:, 3] = t
        pos_x += linx(aa) * dt
        pos_y += liny(aa) * dt
        pos_z += linz(aa) * dt

    pos_final = np.empty([nx, ny, nz, 3])
    ind = 0
    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):
                pos_final[i, j, k, :] = pos_x[ind], pos_y[ind], pos_z[ind]
                ind += 1
    plot = 0
    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax1.imshow(pos_final[:, :, int(nz / 2), 0])
        ax2 = fig.add_subplot(312)
        ax2.imshow(pos_final[:, :, int(nz / 2), 1])
        ax3 = fig.add_subplot(313)
        ax3.imshow(pos_final[:, :, int(nz / 2), 2])
        plt.show()

    print 'advection computed with brio'
    print 'let\' compute the muthafucking CG tenseor now'

    # interpolator for the final positions (aka flow map)
    a = np.array([minx, miny, minz])  # lower boundaries
    b = np.array([maxx, maxy, maxz])  # upper boundaries

    orders = np.array([nx, ny, nz])
    if cubicCG:
        FMx = CubicSpline(a, b, orders, pos_final[:, :, :, 0])
        FMy = CubicSpline(a, b, orders, pos_final[:, :, :, 1])
        FMz = CubicSpline(a, b, orders, pos_final[:, :, :, 2])
    else:
        FMx = LinearSpline(a, b, orders, pos_final[:, :, :, 0])
        FMy = LinearSpline(a, b, orders, pos_final[:, :, :, 1])
        FMz = LinearSpline(a, b, orders, pos_final[:, :, :, 2])

        pos_x, pos_y, pos_z = np.mgrid[minx:maxx - dx:nx * 1j, miny:maxy - dx:ny * 1j,
                              minz:maxz - dx:nz * 1j]

    pos_x = pos_x.ravel()
    pos_y = pos_y.ravel()
    pos_z = pos_z.ravel()
    delta = 1e-4

    ap = np.empty([pos_x.shape[0], 3])
    ap[:, 0] = pos_x + delta
    ap[:, 1] = pos_y
    ap[:, 2] = pos_z
    am = np.empty([pos_x.shape[0], 3])
    am[:, 0] = pos_x - delta
    am[:, 1] = pos_y
    am[:, 2] = pos_z
    g00 = (FMx(ap) - FMx(am)) / (2 * delta)
    g10 = (FMy(ap) - FMy(am)) / (2 * delta)
    g20 = (FMz(ap) - FMz(am)) / (2 * delta)

    ap = np.empty([pos_x.shape[0], 3])
    ap[:, 0] = pos_x
    ap[:, 1] = pos_y + delta
    ap[:, 2] = pos_z
    am = np.empty([pos_x.shape[0], 3])
    am[:, 0] = pos_x
    am[:, 1] = pos_y - delta
    am[:, 2] = pos_z
    g01 = (FMx(ap) - FMx(am)) / (2 * delta)
    g11 = (FMy(ap) - FMy(am)) / (2 * delta)
    g21 = (FMz(ap) - FMz(am)) / (2 * delta)

    ap = np.empty([pos_x.shape[0], 3])
    ap[:, 0] = pos_x
    ap[:, 1] = pos_y
    ap[:, 2] = pos_z + delta
    am = np.empty([pos_x.shape[0], 3])
    am[:, 0] = pos_x
    am[:, 1] = pos_y
    am[:, 2] = pos_z - delta
    g02 = (FMx(ap) - FMx(am)) / (2 * delta)
    g12 = (FMy(ap) - FMy(am)) / (2 * delta)
    g22 = (FMz(ap) - FMz(am)) / (2 * delta)

    ggrid = np.empty([nx, ny, nz, 3, 3])
    ind = 0
    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):
                ggrid[i, j, k, 0, 0] = g00[ind]
                ggrid[i, j, k, 0, 1] = g01[ind]
                ggrid[i, j, k, 0, 2] = g02[ind]
                ggrid[i, j, k, 1, 0] = g10[ind]
                ggrid[i, j, k, 1, 1] = g11[ind]
                ggrid[i, j, k, 1, 2] = g12[ind]
                ggrid[i, j, k, 2, 0] = g20[ind]
                ggrid[i, j, k, 2, 1] = g21[ind]
                ggrid[i, j, k, 2, 2] = g22[ind]
                ind += 1
    plot = 0
    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax1.imshow(ggrid[:, :, int(nz / 2), 0, 0])
        ax2 = fig.add_subplot(312)
        ax2.imshow(ggrid[:, :, int(nz / 2), 1, 1])
        ax3 = fig.add_subplot(313)
        ax3.imshow(ggrid[:, :, int(nz / 2), 2, 2])
        plt.show()

    print " CG tensor computed"
    eigenValues, eigenVectors = LA.eig(ggrid)

    eigvec1 = np.empty((nx, ny, nz, 3))
    eigvec3 = np.empty((nx, ny, nz, 3))
    eigval1 = np.empty((nx, ny, nz))
    eigval3 = np.empty((nx, ny, nz))
    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):
                mineig = np.argmin(eigenValues[i, j, k, :])
                maxeig = np.argmax(eigenValues[i, j, k, :])
                eigval1[i, j, k] = eigenValues[i, j, k, mineig]
                eigval3[i, j, k] = eigenValues[i, j, k, maxeig]
                eigvec1[i, j, k, :] = eigenVectors[i, j, k, :, mineig]
                eigvec3[i, j, k, :] = eigenVectors[i, j, k, :, maxeig]
    print " eigVal, eigVec computed"
    plot = 0
    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.imshow(eigval3[:, :, int(nz / 2)])
        ax2 = fig.add_subplot(212)
        ax2.imshow(eigval1[:, :, int(nz / 2)])
        plt.show()

    return ggrid, eigval1, eigval3, eigvec1, eigvec3
