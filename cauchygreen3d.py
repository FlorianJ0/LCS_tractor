__author__ = 'p0054421'
__date__ = '30 Juillet 2015'

import numpy as np
from scipy import *
# from vtk.util import numpy_support as vn
# import pygsl._numobj as numx
import time
from pygsl import odeiv
# from pygsl import spline, errors
# from pygsl import _numobj as numx
import scipy.version as ver
from numpy import linalg as LA

import matplotlib.pyplot as plt
from Scientific.Functions.Interpolation import InterpolatingFunction as IF
# from scipy.interpolate import RegularGridInterpolator as IF
from scipy.integrate import ode
from scipy.interpolate import griddata
import inpolator
from scipy.ndimage.filters import gaussian_filter
import sys
from airkaeffe import rk45, heun, euler
sys.path.append('pytricubic-master/')
import tricubic


# imp.load_source('pytricubic-master/tricubic.so')
# import vtk
# from vtk import *

# noinspection PyPep8Naming

def rpog():
    print("# "),


def cgstki3(velp, zplan, tt, dt, nx, ny, nz, dim, domain, simtstep):
    # x = np.arange(tt)
    ttt = int(tt / dt)
    ttt = np.arange(ttt)
    # print ttt[0], ttt[-1]
    n = len(ttt)
    N = 25  # rk45 int step
    integrator = 'dopri5'  # ou dopri5 pour du dormant-prince rk45a-
    rrk45 = 2  # 0=rk45, 1=heun, 2= euler

    # small (so is your dick) vector d1(d1 0 0) d2(0 d2 0) d3(0 0 d3)
    d1 = 0.9
    d2 = 0.9
    d3 = 0.9

    # tranche = zplan  # index de la tranche evaluee
    integ = 'rk45'

    # if doublegyre
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print 'WARING DOMAIN CUT FOR TEST PURPOSE'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    # velp = velp[:, :60, :, :, :]
    # ny = velp.shape[1]
    #
    # print 'velp shape'
    # print velp.shape

    # ptlist = np.indices((nx, ny, nz))
    # ptlist = ptlist.astype(float)
    domain = domain.astype(float)
    dx = abs(domain[1] - domain[0]) / nx
    dy = abs(domain[3] - domain[2]) / ny
    dz = abs(domain[5] - domain[4]) / ny
    print 'dx', dx
    print 'dy', dy
    print 'dz', dz
    rr = 1

    nnx = nx * rr
    nny = ny * rr
    nnz = nz * rr  # pas de sous divisions sur z
    ddx = dx / rr
    ddy = dy / rr
    ddz = dz / rr
    zzplan = zplan * rr
    print "domain:", domain
    print 'zplan =', zplan  # * ddz + domain[4]
    tranche = int((zplan - domain[4]) / ((domain[5] - domain[4]) / nz))

    print 'tranche evaluee %i' % tranche

    print 'dt', dt, 't physique', tt, '# time steps', ttt

    # interp spatiale sur une grille r fois plus fine


    grid = np.indices((nnx, nny, nnz))
    grid = grid.astype(float)

    grid_i = np.empty((dim, nnx, nny, nnz))
    grid_i[0, :, :] = grid[0]
    grid_i[1, :, :] = grid[1]
    grid_i[2, :, :] = grid[2]
    # print grid_i
    # grid_iini = np.zeros((dim, nnx, nny, nnz))
    grid_iini = np.empty((dim, nnx, nny, nnz))
    grid_iini[0, :, :] = grid[0]
    grid_iini[1, :, :] = grid[1]
    grid_iini[2, :, :] = grid[2]

    print 'grind ini shape'
    print grid_iini.shape

    print 'interpolation over space, dx/ %i' % rr
    if np.int(ver.full_version[2:4] < 14):
        print 'scipy version %s 2 low. plz upgrade to 0.14.xx' % ver.full_version
        quit()
    else:
        print 'scipy version %s is high, so am I' % ver.full_version

    # x = np.linspace(domain[0], domain[1], nx)
    # y = np.linspace(domain[2], domain[3], ny)
    # z = np.linspace(domain[4], domain[5], nz)
    x = np.arange(0, nnx, 1)
    y = np.arange(0, nny, 1)
    z = np.arange(0, nnz, 1)

    velpu = velp[:, :, :, 0, :]
    velpv = velp[:, :, :, 1, :]
    velpw = velp[:, :, :, 2, :]
    interpU_i = velp
    stamp = time.time()

    print 'linear interp from Scientific python'
    axes = (x, y, z, ttt)
    fu = IF(axes, velpu)
    fv = IF(axes, velpv)
    fw = IF(axes, velpw)

    print 'avection time !'

    # eq dif => x point = v(x,t)
    # ou: y point = f(x,t)
    # fonction f:
    def f_u(yy, t):
        a = np.array(
            [fu(yy[0], yy[1], yy[2], t) / ddx, fv(yy[0], yy[1], yy[2], t) / ddy, fw(yy[0], yy[1], yy[2], t) / ddz])
        return a

    # print 'fu'. f_u()
    solver = ode(f_u)
    solver.set_integrator('dopri5', rtol=0.001, atol=1e-3)
    t = np.linspace(0, tt, N)
    if rrk45 == 0:
        print 'rk45, a checker'
        toto = 0
        quit()
        for i in xrange(nnx):
            for j in xrange(nny):
                y0 = grid_iini[:, i, j, tranche]
                # if np.all(np.abs(f_u(grid_i[:, i, j, tranche],0)) > np.array([1e-7,1e-7,1e-7])):
                if np.all(np.abs(velp[i, j, tranche, :]) > np.array([1e-7, 1e-7, 1e-7])):
                    grid_i[:, i, j, tranche], err = rk45(f_u, y0, t)[-1]
                    if np.max(err) > 1e-5:
                        print err, 'erreur d integ trop grande, my nigga'
                else:
                    grid_i[:, i, j, k] = [0, 0, 0]  # grid_i[:, i, j, tranche] = y0
                    # on en profite pour faire le mask!!!

                    toto += 1
        print '%i point skipped, ie. %f percents of total points of the slice' % (toto, 100 * toto / (ny * nx))
    elif rrk45 == 1:
        print 'heun'
        toto = 0
        for i in xrange(nnx):
            for j in xrange(nny):
                for k in range(tranche - 2, tranche + 3):
                    y0 = grid_iini[:, i, j, k]
                    if np.all(np.abs(velp[i, j, k, :, 0]) > np.array([1e-7, 1e-7, 1e-7])):
                        grid_i[:, i, j, k] = heun(f_u, y0, t)[-1]
                    else:
                        grid_i[:, i, j, k] = [0, 0, 0]
                        toto += 1
        print '%i point skipped, ie. %f percents of total points of the slice' % (toto, 100 * toto / (ny * nx))
    elif rrk45 == 2:
        print 'euler'
        toto = a = 0
        print ('0            50          100%')
        for i in xrange(nnx):
            for j in xrange(nny):
                for k in range(tranche - 1, tranche + 2):
                    y0 = grid_iini[:, i, j, k]

                    if np.all(np.abs(velp[i, j, k, :, 0]) > np.array([1e-7, 1e-7, 1e-7])):
                        grid_i[:, i, j, k] = euler(f_u, y0, t)[-1]
                    else:
                        # grid_i[:, i, j, k] = [0, 0, 0]
                        toto += 1
            a = 1. * i / nnx
            if (100 * a) % 10 < 1e-3:
                rpog()

        print '\n %i point skipped, ie. %f percents of total points of the domain' % (toto, 100 * toto / (ny * nx * 5))
    else:
        print 'wut ?'
        quit()

    # else:
    #     print grid_i[:,35,24,37]
    #     print 'totototot'
    #     t = ttt
    #     t0 = t[0]
    #     t1 = 0.001#t[-1]
    #     i = j = 0
    #     for i in xrange(nnx):
    #         for j in xrange(nny):
    #             y0 = grid_i[:, i, j, tranche]
    #             solver.set_initial_value(y0, t0)
    #             sol = np.empty((N, 3))
    #             sol[0] = y0
    #             k = 1
    #             dt=1e-3
    #             while solver.successful() and solver.t < t1:
    #                 solver.integrate(t[k])
    #                 sol[k] = solver.y
    #                 k += 1
    #             grid_i[:, i, j, tranche] = sol[-1, :]


    print '-----------------------------------------------------'
    print 'Velocity advected  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'

    FTF = True
    if FTF:
        stamp = time.time()
        # gradient of the flow map
        # shadden method
        # (u, v, w) sur (x, y)
        dphi = np.empty((nnx, nny, 3, 3))

        tricu = True

        if tricu:
            print 'tricubic interp'
            du = tricubic.tricubic(list(grid_i[0, :, :, :]),
                                   [nnx, nny, nnz])  # initialize interpolator with input data on cubic grid
            dv = tricubic.tricubic(list(grid_i[1, :, :, :]),
                                   [nnx, nny, nnz])  # initialize interpolator with input data on cubic grid
            dw = tricubic.tricubic(list(grid_i[2, :, :, :]),
                                   [nnx, nny, nnz])  # initialize interpolator with input data on cubic grid

            tata = 0
            print ('0            50          100%')

            # 3d version haller ann. rev. fluid 2015
            for i in range(1, nnx - 1):
                for j in range(1, nny - 1):
                    if np.all(np.abs(velp[i, j, k, :, 0]) > np.array([1e-7, 1e-7, 1e-7])):

                        dphi[i, j, 0, 0] = (du.ip(list(np.array([i + d1, j, tranche]))) - du.ip(
                            list(np.array([i - d1, j, tranche])))) / (2 * d1)

                        dphi[i, j, 0, 1] = (du.ip(list(np.array([i, j + d2, tranche]))) - du.ip(
                            list(np.array([i, j - d2, tranche])))) / (2 * d2)
                        dphi[i, j, 0, 2] = (du.ip(list(np.array([i, j, tranche + d3]))) - du.ip(
                            list(np.array([i, j, tranche - d3])))) / (2 * d3)

                        dphi[i, j, 1, 0] = (dv.ip(list(np.array([i + d1, j, tranche]))) - dv.ip(
                            list(np.array([i - d1, j, tranche])))) / (2 * d1)
                        dphi[i, j, 1, 1] = (dv.ip(list(np.array([i, j + d2, tranche]))) - dv.ip(
                            list(np.array([i, j - d2, tranche])))) / (2 * d2)
                        dphi[i, j, 1, 2] = (dv.ip(list(np.array([i, j, tranche + d3]))) - dv.ip(
                            list(np.array([i, j, tranche - d3])))) / (2 * d3)

                        dphi[i, j, 2, 0] = (dw.ip(list(np.array([i + d1, j, tranche]))) - dw.ip(
                            list(np.array([i - d1, j, tranche])))) / (2 * d1)
                        dphi[i, j, 2, 1] = (dw.ip(list(np.array([i, j + d2, tranche]))) - dw.ip(
                            list(np.array([i, j - d2, tranche])))) / (2 * d2)
                        dphi[i, j, 2, 2] = (dw.ip(list(np.array([i, j, tranche + d3]))) - dw.ip(
                            list(np.array([i, j, tranche - d3])))) / (2 * d3)

                    else:
                        dphi[i, j, :, :] = np.zeros((3, 3))
                        tata += 1
                a = 1. * i / nnx
                if (100 * a) % 10 < 1e-3:
                    rpog()
            print '\n', tata, ' skipped of,', nnx*nny

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

            gdphi = np.empty((nnx, nny, 3, 3))
            for i in xrange(nnx):
                for j in xrange(nny):
                    gdphi[i, j, :, :] = np.dot(dphi[i, j, :, :].T, dphi[i, j, :, :])


        else:
            quit()
            # axes = (xx, yy, zz)
            # du = IF(axes, dispu)
            # dv = IF(axes, dispv)
            # dw = IF(axes, dispw)
            #
            # d1 = dx / 3
            # d2 = dy / 3
            # d3 = dz / 3
            #
            # # 3d version haller ann. rev. fluid 2015
            # for i in range(1, nnx - 1):
            #     for j in range(1, nny - 1):
            #         # for k in range(1, nnz - 1):
            #         # ACHTUNG CALCUL 3D MAIS SEED 2D
            #         ii = i * ddx + domain[0]
            #         jj = j * ddy + domain[2]
            #         zzzplan = zplan * dz + domain[4]
            #         # print ii,jj,zzzplan
            #         # print
            #         # print ii, jj, zzzplan
            #
            #         dphi[i, j, 0, 0] = (du(ii + d1, jj, zzzplan) - du(ii - d1, jj, zzzplan)) / (2 * d1)
            #         print dphi[i, j, 0, 0]
            #         dphi[i, j, 0, 1] = (du(ii, jj + d2, zzzplan) - du(ii, jj - d2, zzzplan)) / (2 * d2)
            #         dphi[i, j, 0, 2] = (du(ii, jj, zzzplan + d3) - du(ii, jj, zzzplan - d3)) / (2 * d3)
            #
            #         dphi[i, j, 1, 0] = (dv(ii + d1, jj, zzzplan) - dv(ii - d1, jj, zzzplan)) / (2 * d1)
            #         dphi[i, j, 1, 1] = (dv(ii, jj + d2, zzzplan) - dv(ii, jj - d2, zzzplan)) / (2 * d2)
            #         dphi[i, j, 1, 2] = (dv(ii, jj, zzzplan + d3) - dv(ii, jj, zzzplan - d3)) / (2 * d3)
            #
            #         dphi[i, j, 2, 0] = (dw(ii + d1, jj, zzzplan) - dw(ii - d1, jj, zzzplan)) / (2 * d1)
            #         dphi[i, j, 2, 1] = (dw(ii, jj + d2, zzzplan) - dw(ii, jj - d2, zzzplan)) / (2 * d2)
            #         dphi[i, j, 2, 2] = (dw(ii, jj, zzzplan + d3) - dw(ii, jj, zzzplan - d3)) / (2 * d3)
            #
            # print 'dphi[50,50,1,1]'
            # print dphi[50, 50, 1, 1]
            # # bords a l arrache;
            # dphi[0, :, 0, 0] = dphi[1, :, 0, 0]
            # dphi[nnx - 1, :, 0, 0] = dphi[nnx - 2, :, 0, 0]
            # dphi[:, 0, 0, 0] = dphi[:, 1, 0, 0]
            # dphi[:, nny - 1, 0, 0] = dphi[:, nny - 2, 0, 0]
            #
            # dphi[0, :, 0, 1] = dphi[1, :, 0, 1]
            # dphi[nnx - 1, :, 0, 1] = dphi[nnx - 2, :, 0, 1]
            # dphi[:, 0, 0, 1] = dphi[:, 1, 0, 0]
            # dphi[:, nny - 1, 0, 1] = dphi[:, nny - 2, 0, 1]
            #
            # dphi[0, :, 1, 0] = dphi[1, :, 1, 0]
            # dphi[nnx - 1, :, 1, 0] = dphi[nnx - 2, :, 1, 0]
            # dphi[:, 0, 1, 0] = dphi[:, 1, 1, 0]
            # dphi[:, nny - 1, 1, 0] = dphi[:, nny - 2, 1, 0]
            #
            # dphi[0, :, 1, 1] = dphi[1, :, 1, 1]
            # dphi[nnx - 1, :, 1, 1] = dphi[nnx - 2, :, 1, 1]
            # dphi[:, 0, 1, 1] = dphi[:, 1, 1, 1]
            # dphi[:, nny - 1, 1, 1] = dphi[:, nny - 2, 1, 1]
            #
            # gdphi = np.empty((nnx, nny, 3, 3))
            # for i in xrange(nnx):
            #     for j in xrange(nny):
            #         gdphi[i, j, :, :] = np.dot(dphi[i, j, :, :].T, dphi[i, j, :, :])

        # print dphi.shape, dphi.T.shape,gdphi.shape
        # toto=np.inner(dphi.T, dphi)
        # print np.array_equal(gdphi,toto)
        # print '------------------------'

        eigenValues, eigenVectors = LA.eig(gdphi)

        eigvec1 = np.empty((nnx, nny, 3))
        eigvec3 = np.empty((nnx, nny, 3))
        eigval1 = np.empty((nnx, nny))
        eigval3 = np.empty((nnx, nny))

        # print 'min,max,avg,stdev'
        # a = eigenValues[:, :, 0] * eigenValues[:, :, 1] * eigenValues[:, :, 2]
        # # print np.min(a)
        # # print np.max(a)
        # # print np.average(a)
        # # print np.std(a)
        # # print '{{{{{{'

        for i in xrange(eigenValues.shape[0]):
            for j in xrange(eigenValues.shape[1]):
                eigval1[i, j] = np.min(eigenValues[i, j, :])
                eigval3[i, j] = np.max(eigenValues[i, j, :])
                eigvec1[i, j, :] = eigenVectors[i, j, :, np.argmin(eigenValues[i, j, :])]
                eigvec3[i, j, :] = eigenVectors[i, j, :, np.argmax(eigenValues[i, j, :])]


        print '-----------------------------------------------------'
        print 'Flow map and eigval/eigvec computed in %f s ' % (time.time() - stamp)
        print '-----------------------------------------------------'

    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex=True, sharey=True)
    # print didx.shape
    Y, X = np.mgrid[0:nx * dx:rr * nx * 1j, 0:ny * dy:rr * ny * 1j]

    uu = grid_i[0, :, :, zzplan] - grid_iini[0, :, :, zzplan]  # -grid_iini[0,:,:]
    vv = grid_i[1, :, :, zzplan] - grid_iini[1, :, :, zzplan]  # -grid_iini[1,:,:]
    ww = grid_i[2, :, :, zzplan] - grid_iini[2, :, :, zzplan]  # -grid_iini[1,:,:]
    magx = np.sqrt(uu * uu + vv * vv + ww * ww)
    U = interpU_i[:, :, 0, 0]
    V = interpU_i[:, :, 1, 0]
    magu = np.sqrt(U * U + V * V)
    # print grid_i[0, 5, :]- grid_iini[0, 5, :]
    ax4.imshow(velpu[:, :, tranche, 0], vmin=-0.05, vmax=0.05, cmap='jet', aspect='auto')
    ax5.imshow(velpv[:, :, tranche, 0], vmin=-0.05, vmax=0.05, cmap='jet', aspect='auto')
    ax6.imshow(velpw[:, :, tranche, 0], vmin=-0.05, vmax=0.05, cmap='jet', aspect='auto')
    # ax2.imshow(dispu[:,:,tranche-1]-dispu[:,:,tranche+1])
    # ax3.imshow(dispu[:,:,tranche+1]-grid_iini[0,:,:,tranche+1])
    # ax3.imshow(grid_i[2, :, :,zzplan])
    # ax2.imshow(magx)
    ax1.imshow(grid_i[0, :, :, tranche])
    ax2.imshow(grid_i[1, :, :, tranche])
    ax3.imshow(grid_i[2, :, :, tranche])

    ax7.imshow(dphi[:, :, 0, 0])
    with file('test.txt', 'w') as outfile:
        np.savetxt(outfile, dphi[:, :, 0, 0])
    ax8.imshow(dphi[:, :, 0, 1])
    ax9.imshow(dphi[:, :, 1, 1])
    # ax2.imshow(didy)



    # ax3.quiver(X, Y, U, V, color=magu)
    # ax4.streamplot(X, Y, uu, vv, density=0.6, color='k', linewidth=magx)

    plt.show()

    print '-------------------------'
    print 'error', np.random.random_integers(0, 100)
    print '-------------------------'
    return eigval1, eigval3, eigvec1, eigvec3, interpU_i
    # return
