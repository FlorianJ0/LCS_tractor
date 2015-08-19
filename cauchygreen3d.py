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
import tricubic

# import vtk
# from vtk import *

# noinspection PyPep8Naming
def cgstki3(vel, zplan, tt, dt, nx, ny, nz, dim, domain, simtstep):
    x = np.arange(tt)
    ttt = int(tt / dt)
    ttt = np.arange(ttt)

    n = len(ttt)
    N = 15   #rk45 int step
    integrator = 'vode'#ou dopri5 pour du dormant-prince rk45a-

    tranche = zplan  # index de la tranche evaluee
    velp = np.empty((nx, ny, nz, dim, n))
    # velp[x,y,(i,j,k), t]

    # velp = vel

    # on recup une tranche du domaine et les 2 composantes de vitesse associees
    integ = 'rk45'

    # if doublegyre
    print 'vel shape'
    print vel.shape
    velp = vel
    # velp= gaussian_filter(velp, 0.5)
    # print velp[25,25,:,:]

    ptlist = np.indices((nx, ny, nz))
    ptlist = ptlist.astype(float)
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
    print 'zplan =', zplan * ddz + domain[4]

    ptlist[0] = ptlist[0] * dx + domain[0]
    ptlist[1] = ptlist[1] * dy + domain[2]
    ptlist[2] = ptlist[2] * dz + domain[4]
    # ptlistini = ptlist

    # stamp = time.time()
    # interpolator
    h = dt  # lenght t step
    # n = int((tt) / h)  # nb t step
    # n=tt
    # n += 1
    # print tt, n
    # print '%i interpolation onver time intervals' % n
    # noinspection PyPep8Naming

    print 'dt', dt, 't physique', tt, '# time steps', ttt

    # interp spatiale sur une grille r fois plus fine


    grid = np.indices((nnx, nny, nnz))
    grid = grid.astype(float)
    # grid=np.swapaxes(grid,1,2)
    grid[0] = grid[0] * abs(domain[1] - domain[0]) / nnx + domain[0]
    grid[1] = grid[1] * abs(domain[3] - domain[2]) / nny + domain[2]
    grid[2] = grid[2] * abs(domain[5] - domain[4]) / nnz + domain[4]
    # gridini = grid



    # interpU_i = np.empty((nnx, nny, nnz, dim, n))
    # newgridu = np.empty((nnx, nny, nnz, n))
    # newgridv = np.empty((nnx, nny, nnz, n))
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

    # print grid[0].shape
    # print interpU[25, 25, 0, 0]
    # print grid_i
    # timeinter = False
    # if timeinter:
    # a putain de refaire avec inperpolate.interpnd
    # k = 0
    print 'interpolation over space, dx/ %i' % rr
    if np.int(ver.full_version[2:4] < 14):
        print 'scipy version %s 2 low. plz upgrade to 0.14.xx' % ver.full_version
        quit()
    else:
        print 'scipy version %s is high, so am I' % ver.full_version


    # inpolator.interp(nx,ny,nz,velp,0,0, rr)
    # print velp[:,:,:,0,0].shape
    # quit()
    # interp tranche par tranche, parceque c'est ce qui  nous interesse dumbass
    # print domain
    # print nx, ny, nz, nnx, nny, nnz
    # points = np.empty((nx * ny, 3))
    # val = np.empty((nx * ny, 3, n))
    x = np.linspace(domain[0], domain[1], nx)
    y = np.linspace(domain[2], domain[3], ny)
    z = np.linspace(domain[4], domain[5], nz)
    velpu = velp[:, :, :, 0, :]
    velpv = velp[:, :, :, 1, :]
    velpw = velp[:, :, :, 2, :]





    # print pts
    # for t in ttt:
        # vals = np.array([velpu[:,0,0,0,t],velp[:,0,0,0,t],velp[1]]).swapaxes(1,0)
    # t = ttt
    # print z, 'nfABDASdbashjkdb'
    # print 'nz', nz, domain[4],domain[5]
    # x_i = np.linspace(domain[0],domain[1], nnx)
    # y_i = np.linspace(domain[2],domain[3], nny)
    # z_i = np.linspace(domain[4],domain[5], nnz)



    axes = (x, y, z, ttt)
    fu = IF(axes, velpu)
    fv = IF(axes, velpv)
    fw = IF(axes, velpw)

    interpU_i = velp

    print 'avection time !'
    # print grid_i
    euler=False
    rk45=True
    if euler:
        print grid_i.shape
        print grid_i[0, 0, 0, zplan]
        print 'ptiny'
        for ti in ttt:
            print 'advection from time ', ti * dt, 'to ', dt * (ti + 1)
            print 'ti', ti
            for i in xrange(nnx):
                for j in xrange(nny):
                    xxx = np.array([grid_i[0, i, j, zzplan], grid_i[1, i, j, zzplan], grid_i[2, i, j, zzplan]])
                    grid_i[0, i, j, zzplan] += fu(xxx[0], xxx[1], xxx[2], ti * dt) * dt
                    grid_i[1, i, j, zzplan] += fv(xxx[0], xxx[1], xxx[2], ti * dt) * dt
                    grid_i[2, i, j, zzplan] += fw(xxx[0], xxx[1], xxx[2], ti * dt) * dt
            # print grid_i[0,50,50,25]*grid_i[0,50,50,25]+grid_i[1,50,50,25]*grid_i[1,50,50,25]+grid_i[2,50,50,25]*grid_i[2,50,50,25]

        print grid_i.shape
        print 'grid_i.shape'
    elif rk45:
        #eq dif => x point = v(x,t)
        # ou: y point = f(x,t)
        #fonction f:
        def f_u(t, y):
            return np.array([fu(y[0],y[1],y[2], t),fv(y[0],y[1],y[2], t),fw(y[0],y[1],y[2], t)])

        #
        # print velp[25,25,20,0,0]
        # print fu(25*ddx,25*ddy,20*ddz,1*dt)

        # r = odeint(f_u).set_integrator('dopri5')
        solver = ode(f_u)
        solver.set_integrator(integrator)
        # t1=ttt[-1]*dt
        # y0=np.array([0.1,0.1,0.25])
        # print f_u(0,y0)
        print '-----------------------------'
        # t0=0
        # y0=grid_iini[:,50,50,25]
        # print fu(y0[0],y0[1],y0[2], 0.2)
        # print fv(y0[0],y0[1],y0[2], 0.2)
        # print fw(y0[0],y0[1],y0[2], 0.2)
        print 'totototot'
        stamp = time.time()
        t = np.linspace(ttt[0]*dt, (ttt[-1])*dt, N)
        t0 = t[0]
        t1 = t[-1]

        for i in xrange(nnx):
            for j in xrange(nny):
                y0=grid_iini[:,i,j,zzplan]
                solver.set_initial_value(y0, t0)
                # print y0
                sol = np.empty((N, 3))
                sol[0] = y0
                k = 1
                while solver.successful() and solver.t < t1:
                    solver.integrate(t[k])
                    sol[k] = solver.y
                    k += 1
                grid_i[:, i, j, zzplan] = sol[-1,:]
                # print sol[0,:]
                # print sol[-1,:]
                # print 'tototo'
                # print grid_i[:, i, j, zzplan]
        print 'advection time %f s ' %(time.time()-stamp)

        # print grid_i[:, :, :, zzplan]
        # print 'y0', y0
        # print 't0', t0
        # print 't',t


        # print sol[:,0]*sol[:,0]+sol[:,1]*sol[:,1]+sol[:,2]*sol[:,2]
        # r.set_initial_value(y0, t0).set_f_params(0)
        # while r.successful() and r.t < t1:
        #     r.integrate(r.t+dt)
        #     print("%g %g" % (r.t, r.y))
        print 'rk45'


    FTF = True
    if FTF:
        # gradient of the flow map
        # shadden method
        # (u, v, w) sur (x, y)
        dphi = np.empty((nnx, nny, 3, 3))
        # 0,0 0,1 0,2
        # 1,0 1,1 1,2
        # 2,0 2,1 2,2
        xx = np.linspace(domain[0], domain[1], nnx)
        yy = np.linspace(domain[2], domain[3], nny)
        zz = np.linspace(domain[4], domain[5], nnz)
        #
        dispu = grid_i[0, :, :, :]
        dispv = grid_i[1, :, :, :]
        dispw = grid_i[2, :, :, :]
        #
        axes = (xx, yy, zz)
        du = IF(axes, dispu)
        dv = IF(axes, dispv)
        dw = IF(axes, dispw)

        d1 = dx / 3
        d2 = dy / 3
        d3 = dz / 3

        # 3d version haller ann. rev. fluid 2015
        for i in range(1, nnx - 1):
            for j in range(1, nny - 1):
                # for k in range(1, nnz - 1):
                # ACHTUNG CALCUL 3D MAIS SEED 2D
                ii = i * ddx + domain[0]
                jj = j * ddy + domain[2]
                zzzplan = zplan * dz + domain[4]
                # print ii,jj,zzzplan
                # print
                # print ii, jj, zzzplan

                dphi[i, j, 0, 0] = (du(ii + d1, jj, zzzplan) - du(ii - d1, jj, zzzplan)) / (2 * d1)
                dphi[i, j, 0, 1] = (du(ii, jj + d2, zzzplan) - du(ii, jj - d2, zzzplan)) / (2 * d2)
                dphi[i, j, 0, 2] = (du(ii, jj, zzzplan + d3) - du(ii, jj, zzzplan - d3)) / (2 * d3)

                dphi[i, j, 1, 0] = (dv(ii + d1, jj, zzzplan) - dv(ii - d1, jj, zzzplan)) / (2 * d1)
                dphi[i, j, 1, 1] = (dv(ii, jj + d2, zzzplan) - dv(ii, jj - d2, zzzplan)) / (2 * d2)
                dphi[i, j, 1, 2] = (dv(ii, jj, zzzplan + d3) - dv(ii, jj, zzzplan - d3)) / (2 * d3)

                dphi[i, j, 2, 0] = (dw(ii + d1, jj, zzzplan) - dw(ii - d1, jj, zzzplan)) / (2 * d1)
                dphi[i, j, 2, 1] = (dw(ii, jj + d2, zzzplan) - dw(ii, jj - d2, zzzplan)) / (2 * d2)
                dphi[i, j, 2, 2] = (dw(ii, jj, zzzplan + d3) - dw(ii, jj, zzzplan - d3)) / (2 * d3)


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

        # print dphi.shape, dphi.T.shape,gdphi.shape
        # toto=np.inner(dphi.T, dphi)
        # print np.array_equal(gdphi,toto)
        # print '------------------------'

        eigenValues, eigenVectors = LA.eig(gdphi)

        eigvec1 = np.empty((nnx, nny, 3))
        eigvec3 = np.empty((nnx, nny, 3))
        eigval1 = np.empty((nnx, nny))
        eigval3 = np.empty((nnx, nny))

        print 'min,max,avg,stdev'
        a = eigenValues[:,:,0]*eigenValues[:,:,1]*eigenValues[:,:,2]
        print np.min(a)
        print np.max(a)
        print np.average(a)
        print np.std(a)
        print '{{{{{{'

        for i in xrange(eigenValues.shape[0]):
            for j in xrange(eigenValues.shape[1]):
                eigval1[i, j] = np.min(eigenValues[i, j, :])
                eigval3[i, j] = np.max(eigenValues[i, j, :])
                eigvec1[i, j, :] = eigenVectors[i, j, :, np.argmin(eigenValues[i, j, :])]
                eigvec3[i, j, :] = eigenVectors[i, j, :, np.argmax(eigenValues[i, j, :])]


    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    # print didx.shape
    Y, X = np.mgrid[0:nx * dx:rr * nx * 1j, 0:ny * dy:rr * ny * 1j]

    uu = grid_i[0, :, :,zzplan]-grid_iini[0,:,:,zzplan]  # -grid_iini[0,:,:]
    vv = grid_i[1, :, :,zzplan]-grid_iini[1,:,:,zzplan]  # -grid_iini[1,:,:]
    ww = grid_i[2, :, :,zzplan]-grid_iini[2,:,:,zzplan]  # -grid_iini[1,:,:]
    magx = np.sqrt(uu * uu + vv * vv + ww * ww)
    U = interpU_i[:, :, 0, 0]
    V = interpU_i[:, :, 1, 0]
    magu = np.sqrt(U * U + V * V)
    # print grid_i[0, 5, :]- grid_iini[0, 5, :]
    ax1.imshow(velpu[:,:,35,0])
    ax2.imshow(velpv[:,:,35,0])
    ax3.imshow(eigval1)
    # ax2.imshow(magx)
    # ax2.imshow(didy)
    # ax2.imshow(didy)
    ax4.imshow(eigval3)
    # ax3.quiver(X, Y, U, V, color=magu)
    # ax4.streamplot(X, Y, uu, vv, density=0.6, color='k', linewidth=magx)
    plt.show()

    print '-------------------------'
    print 'error', np.random.random_integers(0, 100)
    print '-------------------------'
    return eigval1, eigval3, eigvec1, eigvec3, interpU_i
