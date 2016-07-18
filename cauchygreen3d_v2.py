__author__ = 'p0054421'
__date__ = '30 Juillet 2015'

from numba import jit
from scipy import *

# from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from interpolation.splines import LinearSpline, CubicSpline
import ConfigParser

# noinspection PyPep8Naming
Config = ConfigParser.ConfigParser()
Config.read('parameters.ini')
# cubic interp ? if 0 then linear
cubic = 1


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
    if cubic:
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

    FMx = CubicSpline(a, b, orders, pos_final[:, :, :, 0])
    FMy = CubicSpline(a, b, orders, pos_final[:, :, :, 1])
    FMz = CubicSpline(a, b, orders, pos_final[:, :, :, 2])

    @jit
    def CG():
        ggrid = np.empty([nx, ny, nz, 3, 3])
        delta = 1e-4
        for i in xrange(nx):
            x = minx + i * dx
            for j in xrange(ny):
                y = miny + j * dy
                for k in xrange(nz):
                    z = minz + k * dz
                    ggrid[i, j, k, 0, 0] = FMx(np.array([x + delta, y, z])) - FMx(np.array([x - delta, y, z]))
                    ggrid[i, j, k, 0, 1] = FMx(np.array([x, y + delta, z])) - FMx(np.array([x, y - delta, z]))
                    ggrid[i, j, k, 0, 2] = FMx(np.array([x, y, z + delta])) - FMx(np.array([x, y, z - delta]))

                    ggrid[i, j, k, 1, 0] = FMy(np.array([x + delta, y, z])) - FMy(np.array([x - delta, y, z]))
                    ggrid[i, j, k, 1, 1] = FMy(np.array([x, y + delta, z])) - FMy(np.array([x, y - delta, z]))
                    ggrid[i, j, k, 1, 2] = FMy(np.array([x, y, z + delta])) - FMy(np.array([x, y, z - delta]))

                    ggrid[i, j, k, 2, 0] = FMz(np.array([x + delta, y, z])) - FMz(np.array([x - delta, y, z]))
                    ggrid[i, j, k, 2, 1] = FMz(np.array([x, y + delta, z])) - FMz(np.array([x, y - delta, z]))
                    ggrid[i, j, k, 2, 2] = FMz(np.array([x, y, z + delta])) - FMz(np.array([x, y, z - delta]))

        ggrid /= 2 * delta
        return ggrid

    ggrid = CG()
    plot = 1
    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax1.imshow(ggrid[:, :, int(nz / 2), 0, 0])
        ax2 = fig.add_subplot(312)
        ax2.imshow(ggrid[:, :, int(nz / 2), 1, 1])
        ax3 = fig.add_subplot(313)
        ax3.imshow(ggrid[:, :, int(nz / 2), 2, 2])
        plt.show()

    #


    return
