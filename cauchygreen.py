__author__ = 'p0054421'
__date__ = '30 Juillet 2015'

import numpy as np
from scipy import *
# from vtk.util import numpy_support as vn
# import pygsl._numobj as numx
import time
# from pygsl import odeiv
# from pygsl import spline, errors
# from pygsl import _numobj as numx
from scipy.interpolate import interp1d
from numpy import linalg as LA
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline


# import vtk
# from vtk import *

# noinspection PyPep8Naming
def cgstki(vel, z, tt, dt, nx, ny, dim, domain, simtstep):
    x = np.arange(tt)
    ttt = int(tt / dt)
    ttt = np.arange(ttt)
    n = len(ttt)
    velp = np.empty((nx, ny, dim, n))
    # velp[x,y,(i,j,k), t]

    # velp = vel

    # on recup une tranche du domaine et les 2 composantes de vitesse associees
    integ = 'rk45'

    # if doublegyre
    velp[:, :, :, :] = vel[:, :, z, 0:2, :]
    print velp[25, 25, :, :]

    ptlist = np.indices((nx, ny))
    ptlist = ptlist.astype(float)
    domain = domain.astype(float)
    dx = abs(domain[1] - domain[0]) / nx
    dy = abs(domain[3] - domain[2]) / ny
    print 'dx', dx
    print 'dy', dy
    rr = 1
    nnx = nx * rr
    nny = ny * rr

    ptlist[0] = ptlist[0] * dx + domain[0]
    ptlist[1] = ptlist[1] * dy + domain[2]
    ptlistini = ptlist

    stamp = time.time()
    # interpolator
    h = dt  # lenght t step
    # n = int((tt) / h)  # nb t step
    # n=tt
    # n += 1
    # print tt, n
    # print '%i interpolation onver time intervals' % n
    # noinspection PyPep8Naming

    print 'dt', dt, 't physique', tt, '# time steps', ttt
    interpU = np.empty((nx, ny, dim, n))

    interTime = False
    if interTime:
        for i in xrange(nx):
            for j in xrange(ny):
                # print 'i', i, 'j', j
                y = velp[i, j, :, :]
                # print y.shape
                f = interp1d(x, y, kind='quadratic')
                for ti in ttt:
                    interpU[i, j, :, ti] = f(ti)
    else:
        interpU = velp

    # interp spatiale sur une grille r fois plus fine

    grid = np.indices((nnx, nny))
    grid = grid.astype(float)
    # grid=np.swapaxes(grid,1,2)
    grid[0] = grid[0] * abs(domain[1] - domain[0]) / nnx + domain[0]
    grid[1] = grid[1] * abs(domain[3] - domain[2]) / nny + domain[2]
    gridini = grid

    points = np.empty((nx * ny, dim))
    val = np.empty((nx * ny, dim, n))

    interpU_i = np.empty((nnx, nny, dim, n))
    newgridu = np.empty((nnx, nny, n))
    newgridv = np.empty((nnx, nny, n))
    grid_i = np.empty((dim, nnx, nny))
    grid_i[0, :, :] = grid[0]
    grid_i[1, :, :] = grid[1]
    grid_iini = np.zeros((dim, nnx, nny))
    grid_iini = np.empty((dim, nnx, nny))
    grid_iini[0, :, :] = grid[0]
    grid_iini[1, :, :] = grid[1]
    # print interpU[25, 25, 0, 0]
    # print grid_i
    # timeinter = False
    # if timeinter:
    # a putain de refaire avec inperpolate.interpnd
    # k = 0
    print 'interpolation over space, dx/ %i' % rr
    print domain
    print nx, ny, nnx, nny
    for ti in ttt:
        print 'ti', ti
        k = 0
        for i in xrange(nx):
            for j in xrange(ny):
                points[k] = ptlistini[:, i, j]
                val[k, 0, ti] = interpU[i, j, 0, ti]
                val[k, 1, ti] = interpU[i, j, 1, ti]
                k += 1
        newgridu[:, :, ti] = griddata(points, val[:, 0, ti], (grid[0], grid[1]), method='cubic', fill_value=0.)
        newgridv[:, :, ti] = griddata(points, val[:, 1, ti], (grid[0], grid[1]), method='cubic', fill_value=0.)
        print np.max(newgridv[:, :, ti])
        print '*********************'
    interpU_i[:, :, 0, :] = newgridu[:, :, :]
    # interpU_i[:, :, 0, :] = interpU[:, :, 0,:]
    interpU_i[:, :, 1, :] = newgridv[:, :, :]
    # interpU_i[:, :, 1, :] = interpU[:, :, 1,:]
    print 'toto', points.shape, val.shape, grid.shape
    print np.max(newgridu)
    # print grid[0]
    # else:
    # interpU_i = velp
    print 'avection time !'
    for ti in ttt:
        print 'advection from time ', ti * dt, 'to ', dt * (ti + 1)
        print 'ti', ti
        # totou = interp2d(grid[0, :, :], grid[1, :, :], interpU_i[:, :, 0, ti], kind='linear')
        # print 'x', grid_iini[0, :, 0]
        # print 'y', grid_iini[1, 0, :]
        totou = RectBivariateSpline(grid_iini[0, :, 0], grid_iini[1, 0, :], interpU_i[:, :, 0, ti])
        # totov = interp2d(grid[0, :, :], grid[1, :, :], interpU_i[:, :, 1, ti], kind='linear')
        totov = RectBivariateSpline(grid_iini[0, :, 0], grid_iini[1, 0, :], interpU_i[:, :, 1, ti])
        print 'interpolator interpolated'
        for i in xrange(nnx):
            for j in xrange(nny):
                grid_i[0, i, j] += totou(grid_i[0, i, j], grid_i[1, i, j]) * dt
                grid_i[1, i, j] += totov(grid_i[0, i, j], grid_i[1, i, j]) * dt

    # gradient of the flow map
    # shadden method
    dphi = np.empty((nnx, nny, 2, 2))
    # 0,0 0,1
    # 1,0 1,1
    ftle = True
    if ftle:
        for i in range(1, nnx - 1):
            for j in range(1, nny - 1):
                dphi[i, j, 0, 0] = (grid_i[0, i + 1, j] - grid_i[0, i - 1, j]) / (
                    grid_iini[0, i + 1, j] - grid_iini[0, i - 1, j])
                dphi[i, j, 0, 1] = (grid_i[0, i, j + 1] - grid_i[0, i, j - 1]) / (
                    grid_iini[1, i, j + 1] - grid_iini[1, i, j - 1])
                dphi[i, j, 1, 0] = (grid_i[1, i + 1, j] - grid_i[1, i - 1, j]) / (
                    grid_iini[0, i + 1, j] - grid_iini[0, i - 1, j])
                dphi[i, j, 1, 1] = (grid_i[1, i, j + 1] - grid_i[1, i, j - 1]) / (
                    grid_iini[1, i, j + 1] - grid_iini[1, i, j - 1])

        # bords a l arrache;
        dphi[0, :, 0, 0] = dphi[1, :, 0, 0]
        dphi[nnx - 1, :, 0, 0] = dphi[nnx - 2, :, 0, 0]
        dphi[:, 0, 0, 0] = dphi[:, 1, 0, 0]
        dphi[:, nny - 1, 0, 0] = dphi[:, nny - 2, 0, 0]

        dphi[0, :, 0, 1] = dphi[1, :, 0, 1]
        dphi[nnx - 1, :, 0, 1] = dphi[nnx - 2, :, 0, 1]
        dphi[:, 0, 0, 1] = dphi[:, 1, 0, 0]
        dphi[:, nny - 1, 0, 1] = dphi[:, nny - 2, 0, 1]

        dphi[0, :, 1, 0] = dphi[1, :, 1, 0]
        dphi[nnx - 1, :, 1, 0] = dphi[nnx - 2, :, 1, 0]
        dphi[:, 0, 1, 0] = dphi[:, 1, 1, 0]
        dphi[:, nny - 1, 1, 0] = dphi[:, nny - 2, 1, 0]

        dphi[0, :, 1, 1] = dphi[1, :, 1, 1]
        dphi[nnx - 1, :, 1, 1] = dphi[nnx - 2, :, 1, 1]
        dphi[:, 0, 1, 1] = dphi[:, 1, 1, 1]
        dphi[:, nny - 1, 1, 1] = dphi[:, nny - 2, 1, 1]

        gdphi = np.empty((nnx, nny, 2, 2))
        for i in xrange(nnx):
            for j in xrange(nny):
                gdphi[i, j, :, :] = np.dot(dphi[i, j, :, :].T, dphi[i, j, :, :])
                # gdphi[i, j, :, :] = dphi[i, j, :, :]* dphi[i, j, :, :]

        feuteuleu = np.empty((nnx, nny))
        for i in xrange(nnx):
            for j in xrange(nny):
                feuteuleu[i, j] = np.log(np.sqrt(np.max(LA.eigvals(gdphi[i, j, :, :])[1])))
                pass
                # feuteuleu[i, j] = np.sqrt(LA.eigvals(gdphi[i, j, :, :])[1])
                # print len(LA.eigvals(gdphi[i, j, :, :]))
                # print LA.eigvals(gdphi[20, 20, :,:])

    eigenValues = np.empty((nnx, nny, 2))  # en chaque point, 2 eigval + 2 eigvec
    eigenVectors = np.empty((nnx, nny, 2, 2))  # en chaque point, 2 eigval + 2 eigvec

    print '111'
    eigenValues, eigenVectors = LA.eig(gdphi)

    stamp = time.time()
    for i in xrange(nnx):
        for j in xrange(nny):
            if eigenValues[i, j, 0] < eigenValues[i, j, 0]:
                print i, j, 'EIG NOT ORDERED DAMMIT'
    print 'time= %f' % (time.time() - stamp)

    #
    # f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    # Y, X = np.mgrid[0:nx * dx:rr * nx * 1j, 0:ny * dy:rr * ny * 1j]
    #
    # uu = grid_i[0, :, :] - grid_iini[0, :, :]  # -grid_iini[0,:,:]
    # vv = grid_i[1, :, :] - grid_iini[1, :, :]  # -grid_iini[1,:,:]
    # magx = np.sqrt(uu * uu + vv * vv)
    # U = interpU_i[:, :, 0, 0]
    # V = interpU_i[:, :, 1, 0]
    # magu = np.sqrt(U * U + V * V)
    # ax1.imshow(uu)
    # ax2.imshow(vv)
    # ax3.imshow(feuteuleu)
    # ax4.imshow(magx)


    print '-------------------------'
    print 'error', np.random.random_integers(0, 100)
    print '-------------------------'
    return eigenValues, eigenVectors, interpU_i
