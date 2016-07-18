# python2.7
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

import readdomain
from interpolation.splines import CubicSpline


# def cython_mean(double[:] x):
#     cdef double total = 0
#     for i in range(len(x)):
#         total += x[i]
#     # return total / len(x)
#     return


def green(traj, fold):
    print 'entering gree'
    [nx, ny, nz, dim, tphys, dt, domain] = readdomain.read_files(fold)

    print type(traj)
    print traj.shape
    pos = traj[0]
    vel = traj[1]
    vort = traj[2]
    nPaths = len(vort)
    """
    r = grid refinement coeff for CG tensor computation
    we need intermediary points to compute the CG tenser
    if grid = 0 - 1 - 2 (3 pts, 2 intervals) we need 0 - 0.33 - 0.66 - 1 - 1.33 - 1.66 -2
    (7 pts, 6 intervals)
    npoints for np.mgrid: 3*(n-1)+1
    new pts = nnx...
    """
    r = 1
    dx = (domain[1] - domain[0]) / nx
    dy = (domain[3] - domain[2]) / ny
    dz = (domain[5] - domain[4]) / nz
    if abs(dx - dy) > 1e-4 or abs(dx - dz) > 1e-4:
        print 'bad spacing', dx, dy, dz
        quit()
    start_pts = np.empty([nPaths, 3])
    end_pts = np.empty([nPaths, 3])
    k = 0
    for i in pos[:]:
        start_pts[k] = i[0]
        end_pts[k] = i[-1]
        k += 1

    nnx = r * (nx - 1) + 1
    nny = r * (ny - 1) + 1
    nnz = r * (nz - 1) + 1


    grid_x, grid_y, grid_z = np.mgrid[domain[0]:domain[1] - dx:nnx * 1j, domain[2]:domain[3] - dx:nny * 1j,
                             domain[4]:domain[5] - dx:nnz * 1j]

    ddx = grid_x[1, 0, 0] - grid_x[0, 0, 0]  # new mini dx



    # mappedPos = griddata(start_pts, end_pts, (grid_x, grid_y, grid_z), method='linear')
    mappedPosFile = fold + '/mappedPos'
    print mappedPosFile
    if not os.path.isfile(mappedPosFile + '.npy'):
        print "computing mapped positions"
        t0 = time.time()
        mappedPos = interpolate.griddata(start_pts, end_pts, (grid_x, grid_y, grid_z), method='linear',
                                         fill_value=0)
        np.save(mappedPosFile, mappedPos)
        print 'interpolated in', time.time() - t0, 's'
    else:
        mappedPos = np.load(mappedPosFile + '.npy')
        print "interp file loaded"

    plot = 0
    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax1.imshow(mappedPos[:, :, int(nz / 2), 0])
        ax2 = fig.add_subplot(312)
        ax2.imshow(mappedPos[:, :, int(nz / 2), 1])
        ax3 = fig.add_subplot(313)
        ax3.imshow(mappedPos[:, :, int(nz / 2), 2])
        plt.show()

    """
    run of cgreen tensor computation on original grid (nx ny nz).
    nx = (nnx+2)/3
    we construct the grid of adjacent values, i.e. x+dx, x-dx etx etc
    the new grid will be the shape [[nx,ny,nz],[-dx,+dx,-dy,+dy,-dz,+dz]]
    """
    t0 = time.time()
    ggrid = np.zeros([nx, ny, nz, 3, 3])  # 3x3 matrix at each points
    a = np.array([0,0, 0])  # lower boundaries
    b = np.array([nx, ny, ny])  # upper boundaries
    orders = np.array([nx, ny, nz])  # 10 points along each dimension
    valuesx = mappedPos[:,:,:,0]
    valuesy = mappedPos[:,:,:,1]
    valuesz = mappedPos[:,:,:,2]
    splinex = CubicSpline(a, b, orders, valuesx)
    spliney = CubicSpline(a, b, orders, valuesy)
    splinez = CubicSpline(a, b, orders, valuesz)
    print splinex(np.array([25, 10, 12]))
    print mappedPos[25, 10, 12,0]
    quit()
    delta = 1e-2
    for i in range(1, nx - 1):
        ii = 3 * i
        iip = ii + delta
        iim = ii - delta
        for j in xrange(1, ny - 1):
            jj = 3 * j
            jjp = jj + delta
            jjm = jj - delta
            for k in xrange(1, nz - 1):
                kk = 3 * k
                kkp = kk + delta
                kkm = kk - delta
                ggrid[i, j, k, 0, 0] = splinex([iip, jj, kk]) - splinex([iim, jj, kk])
                ggrid[i, j, k, 0, 1] = splinex([ii, jjp, kk]) - splinex([ii, jjm, kk])
                ggrid[i, j, k, 0, 2] = splinex([ii, jj, kkp]) - splinex([ii, jj, kkm])

                ggrid[i, j, k, 1, 0] = spliney([iip, jj, kk]) - spliney([iim, jj, kk])
                ggrid[i, j, k, 1, 1] = spliney([ii, jjp, kk]) - spliney([ii, jjm, kk])
                ggrid[i, j, k, 1, 2] = spliney([ii, jj, kkp]) - spliney([ii, jj, kkm])

                ggrid[i, j, k, 2, 0] = splinez([iip, jj, kk]) - splinez([iim, jj, kk])
                ggrid[i, j, k, 2, 1] = splinez([ii, jjp, kk]) - splinez([ii, jjm, kk])
                ggrid[i, j, k, 2, 2] = splinez([ii, jj, kkp]) - splinez([ii, jj, kkm])

    ggrid /= 2 * ddx
    print "CG tensor computed in ", time.time() - t0, 's'
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

    return
