__author__ = 'p0054421'
__date__ = '30 Juillet 2015'

from scipy import *
# from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from interpolation.splines import LinearSpline, CubicSpline
import ConfigParser

# noinspection PyPep8Naming
Config = ConfigParser.ConfigParser()
Config.read('parameters.ini')
#cubic interp ? if 0 then linear
cubic = 0


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
    ttt = int(tt / dt)
    ttt = np.arange(ttt)
    tmin = ttt[0] * dt
    tmax = ttt[-1] * dt
    print 'tmin, tmax', tmin, tmax
    # print ttt[0], ttt[-1]
    nt = len(ttt)
    print 'n time step', nt

    # small (so is your dick) vector d1(d1 0 0) d2(0 d2 0) d3(0 0 d3)
    dx = (domain[1] - domain[0]) / nx
    dy = (domain[3] - domain[2]) / ny
    dz = (domain[5] - domain[4]) / nz
    if abs(dx - dy) > 1e-4 or abs(dx - dz) > 1e-4:
        print 'bad spacing', dx, dy, dz
        quit()
    a = np.array([domain[0], domain[2], domain[4], tmin])  # lower boundaries
    b = np.array([domain[1], domain[3], domain[5], tmax])  # upper boundaries

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
    grid_x, grid_y, grid_z, grid_t = np.mgrid[domain[0]:domain[1] - dx:nx * 1j, domain[2]:domain[3] - dx:ny * 1j,
                                     domain[4]:domain[5] - dx:nz * 1j, tmin:tmax:nt * 1j]

    positions = np.vstack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel(), grid_t.ravel()]).T

    # interpUZ = linz(positions)
    # interpUX = linx(positions)
    # interpUY = liny(positions)

    pos_x, pos_y, pos_z = np.mgrid[domain[0]:domain[1] - dx:nx * 1j, domain[2]:domain[3] - dx:ny * 1j,
                          domain[4]:domain[5] - dx:nz * 1j]
    pos_x = pos_x.ravel()
    pos_y = pos_y.ravel()
    pos_z = pos_z.ravel()

    # cheap advection

    t0 = tmin
    z0 = positions
    t1 = tmax
    dt = 0.01  # interp step
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
    plot =
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

    #interpolator for the final positions (aka flow map)
    FMx = CubicSpline(a, b, orders, velpu)
    FMy = CubicSpline(a, b, orders, velpv)
    FMz = CubicSpline(a, b, orders, velpw)

    ggrid = np.empty([nx, ny, nz, 3, 3])
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

    #


    return
'''
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

    # x = np.linspace(0, nnx - 1, nnx, dtype=int16)
    # y = np.linspace(0, nny - 1, nny, dtype=int16)
    # z = np.linspace(0, nnz - 1, nnz, dtype=int16)
    # t = np.linspace(0, len(ttt) - 1, len(ttt), dtype=int16)

    velpu = velp[:, :, :, 0, :]
    velpv = velp[:, :, :, 1, :]
    velpw = velp[:, :, :, 2, :]

    interpU_i = velp
    stamp = time.time()

    # def fu(x, y, z, t):
    #     return velpu[x, y, z, t]
    #
    # def fv(x, y, z, t):
    #     return velpv[x, y, z, t]
    #
    # def fw(x, y, z, t):
    #     return velpw[x, y, z, t]
    #
    # datau = fu(*np.meshgrid(x, y, z, t, indexing='ij', sparse=False))
    # datav = fv(*np.meshgrid(x, y, z, t, indexing='ij', sparse=False))
    # dataw = fw(*np.meshgrid(x, y, z, t, indexing='ij', sparse=False))
    #
    # interpu = RegularGridInterpolator((x, y, z, t), datau, method='nearest', bounds_error=False, fill_value=0)
    # interpv = RegularGridInterpolator((x, y, z, t), datav, method='nearest', bounds_error=False, fill_value=0)
    # interpw = RegularGridInterpolator((x, y, z, t), dataw, method='nearest', bounds_error=False, fill_value=0)


    # quit()
    print 'linear interp NOT from Scientific python'

    print 'avection time !'

    # eq dif => x point = v(x,t)
    # ou: y point = f(x,t)
    # fonction f:
    nn = np.array([nnx, nny, nnz])
    invddx = 1 / ddx
    invddy = 1 / ddy
    invddz = 1 / ddz

    def f_u(yy, t):
        coord = np.array([yy[0], yy[1], yy[2], t])
        # a = np.array(
        #     [interpu(coord) * invddx, interpv(coord) * invddy, interpw(coord) * invddz])

        u = ndimage.map_coordinates(velpu, [[yy[0]], [yy[1]], [yy[2]], [t]], order=3, mode='constant', cval=0.0,
                                    prefilter=False) * invddx
        v = ndimage.map_coordinates(velpv, [[yy[0]], [yy[1]], [yy[2]], [t]], order=3, mode='constant', cval=0.0,
                                    prefilter=False) * invddy
        w = ndimage.map_coordinates(velpw, [[yy[0]], [yy[1]], [yy[2]], [t]], order=3, mode='constant', cval=0.0,
                                    prefilter=False) * invddz
        return np.array([u, v, w])[:, 0]

    # print 'fu'. f_u()
    solver = ode(f_u)
    solver.set_integrator('dopri5', rtol=0.001, atol=1e-3)
    t = np.linspace(0, tt, N)
    bobol = np.zeros((nnx, nny))
    if rrk45 == 0:
        print 'rk45 intregration method, a checker'
        toto = 0
        quit()
        for i in xrange(nnx):
            for j in xrange(nny):
                y0 = grid_iini[:, i, j, tranche]
                # if np.all(np.abs(f_u(grid_i[:, i, j, tranche],0)) > np.array([1e-7,1e-7,1e-7])):
                if np.all(np.abs(velp[i, j, tranche, :]) > np.array([1e-7, 1e-7, 1e-7])):
                    bobol[i, j] = True
                    grid_i[:, i, j, tranche], err = rk45(f_u, y0, t)[-1]
                    if np.max(err) > 1e-5:
                        print err, 'erreur d integ trop grande, my nigga'
                else:
                    grid_i[:, i, j, k] = [0, 0, 0]  # grid_i[:, i, j, tranche] = y0
                    # on en profite pour faire le mask!!!

                    toto += 1
        print '%i point skipped, ie. %f percents of total points of the slice' % (toto, 100 * toto / (ny * nx))
    elif rrk45 == 1:
        print 'heun intregration method'
        toto = 0
        for i in xrange(nnx):
            for j in xrange(nny):
                for k in range(tranche - 2, tranche + 3):
                    y0 = grid_iini[:, i, j, k]
                    if np.all(np.abs(velp[i, j, k, :, 0]) > np.array([1e-7, 1e-7, 1e-7])):
                        bobol[i, j] = True
                        grid_i[:, i, j, k] = heun(f_u, y0, t, nx, ny)[-1]
                    else:
                        grid_i[:, i, j, k] = [0, 0, 0]
                        toto += 1
        print '%i point skipped, ie. %f percents of total points of the slice' % (toto, 100 * toto / (ny * nx))
    elif rrk45 == 2:
        print 'euler intregration method'
        toto = a = 0
        print ('0            50          100%')

        for i in xrange(nnx):
            for j in xrange(nny):
                for k in range(tranche - 0, tranche + 1):
                    y0 = grid_iini[:, i, j, k]

                    if np.all(np.abs(velp[i, j, k, :, 0]) > np.array([1e-7, 1e-7, 1e-7])):
                        bobol[i, j] = True
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
                    # if np.all(np.abs(velp[i, j, k, :, 0]) > np.array([1e-7, 1e-7, 1e-7])):
                    if bobol[i, j]:

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
            print '\n', tata, ' skipped of,', nnx * nny

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
                # print eigvec3[i, j, :]

        # toto = eigval3.astype(short)
        # juliaStacked = np.dstack([toto])
        # x = np.arange(0, nnx)
        # y = np.arange(0, nny)
        # z = np.arange(0, 2)
        # gridToVTK("./julia", x, y, z, cellData = {'julia': juliaStacked})


        print '-----------------------------------------------------'
        print 'Flow map and eigval/eigvec computed in %f s ' % (time.time() - stamp)
        print '-----------------------------------------------------'

        # p1 = win.addPlot(title="Basic array plotting", y=np.random.normal(size=100))














        f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex=True, sharey=True)
        # print didx.shape
        # Y, X = np.mgrid[0:nx * dx:rr * nx * 1j, 0:ny * dy:rr * ny * 1j]

        uu = grid_i[0, :, :, zzplan] - grid_iini[0, :, :, zzplan]  # -grid_iini[0,:,:]
        vv = grid_i[1, :, :, zzplan] - grid_iini[1, :, :, zzplan]  # -grid_iini[1,:,:]
        ww = grid_i[2, :, :, zzplan] - grid_iini[2, :, :, zzplan]  # -grid_iini[1,:,:]
        magx = np.sqrt(uu * uu + vv * vv + ww * ww)
        # U = interpU_i[:, :, 0, 0]
        # V = interpU_i[:, :, 1, 0]
        # magu = np.sqrt(U * U + V * V)
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
        # with file('test.txt', 'w') as outfile:
        # np.savetxt(outfile, dphi[:, :, 0, 0])
        ax8.imshow(dphi[:, :, 0, 1])
        ax9.imshow(dphi[:, :, 1, 1])
        # ax2.imshow(didy)
        # ax3.quiver(X, Y, U, V, color=magu)
        # ax4.streamplot(X, Y, uu, vv, density=0.6, color='k', linewidth=magx)

        plt.show()

    print '-------------------------'
    print 'error', np.random.random_integers(0, 100)
    print '-------------------------'
    return eigval1, eigval3, eigvec1, eigvec3, interpU_i, bobol
    # return
'''