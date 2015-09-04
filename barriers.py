__author__ = 'p0054421'
import numpy as np
import time

import matplotlib.pyplot as plt
from Scientific.Functions.Interpolation import InterpolatingFunction as IF

from airkaeffe import heun


def rpog():
    print("# "),

def helic(vect, dx, dy, dz, eps0):
    print 'vect shape', vect.shape
    nx = vect.shape[0]
    ny = vect.shape[1]
    dd = np.empty([nx, ny, 3])
    ddd = np.empty([nx, ny])
    # eps0 = 5e-1
    # at each point, 3d grad vector of 3D eigenc
    for i in xrange(nx):
        for j in xrange(ny):
            # dd[i, j, :] = np.gradient(vect[i, j, :])
            ddd[i, j] = np.dot(vect[i, j, :], np.gradient(vect[i, j, :]))
            # if abs(ddd[i,j])<eps0:
            #     ddd[i,j]=True
            # else
            #     ddd[i,j]=False

    a = abs(ddd) < eps0
    ddd *= a

    return ddd


def reduced_lines(vect, nx, ny, nz, initpts, thresh, n):
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)
    vecti = (vect[:, :, 0])
    vectj = (vect[:, :, 1])

    # vectj = vect[:,:,1]
    axes = (x, y)
    vi = IF(axes, vecti)
    vj = IF(axes, vectj)

    # vj = IF(axes, vectj)


    def gamma(x, t):
        newpos = np.empty((2))
        # print 'x', x
        a = np.array([vi(x[0], x[1]), vj(x[0], x[1])])
        cond = np.dot(a, np.gradient(a))
        cut = False
        # if (1e-4 < abs(cond) < thresh):
        if (abs(cond) < thresh):
            # i did the maths, must be right
            newpos[0] = - 1. * vj(x[0], x[1])
            newpos[1] = 1. * vi(x[0], x[1])
        else:
            cut = True

        return newpos, cut

    print 'integrate reduced LCSs'
    N = n * 0.2
    # on suppose qu on est toujours normal  a z
    # norm vect = 0 0 -1
    # donc n vectproduct k = kj -ki 0
    # la trajectoire est portee par le vect kj -ki 0 donc dans le plan
    t = np.linspace(0, n, N)  # pseudo integrator
    line = np.zeros((initpts.shape[1], 2, N))
    print ('0            50          100%')
    for i in xrange(initpts.shape[1]):
        y0 = initpts[:, i] * 1.
        line[i, :, :] = heun(gamma, y0, t, nx, ny).swapaxes(1, 0)
        print 'line number %i' % i
        # print y0
        # print line[i, :, -1]
        # print line.shape
        # fname = 'toto' + str(i)
        # np.savetxt(fname, line[i,:,:])
        a = 1. * i / initpts.shape[1]
        if (100 * a) % 10 < 1e-1:
                    rpog()

    print '\n'
    return line


def barrier_type(toto, eigval1, eigval3, eigvec1, eigvec3, vel, tphys, dt, nx, ny, nz, domain, simtstep):
    dx = abs(domain[1] - domain[0]) / nx
    dy = abs(domain[3] - domain[2]) / ny
    dz = abs(domain[5] - domain[4]) / ny
    print 'hyperbolic  LCS'
    # cond helicity
    # H = np.dot(eigvec[0],np.cross(np.gradient() ))
    # eigvecm = eigvec[:,:,0,:]
    # print eigvec.shape
    # print eigvecm.shape
    stamp = time.time()
    initpts3 = helic(eigvec3, dx, dy, dz, 0.005)
    seeds3 = np.array([initpts3.nonzero()[0], initpts3.nonzero()[1]])
    strain_lines = reduced_lines(eigvec3, nx, ny, nz, seeds3, 75,500)
    print '-----------------------------------------------------'
    print 'strain lines (repelling) computed  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'

    stamp = time.time()
    initpts1 = helic(eigvec1, dx, dy, dz, 0.05)
    seeds1 = np.array([initpts1.nonzero()[0], initpts1.nonzero()[1]])
    stretch_lines = reduced_lines(eigvec1, nx, ny, nz, seeds1, 4,5)
    print '-----------------------------------------------------'
    print 'stretch (attracting) lines computed  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'


    ftle = (1/0.16)*np.log(np.sqrt(eigval3+1e-8))
    magU = np.sqrt(vel[:, :, 14, 0, 0] ** 2 + vel[:, :, 14, 1, 0] ** 2 + vel[:, :, 14, 2, 0] ** 2)

    plt.subplot(321)
    plt.imshow(eigval3)
    for j in xrange(strain_lines.shape[0]):
        plt.plot(strain_lines[:, 1, :], strain_lines[:, 0, :], 'w.', ms=1)
    plt.colorbar()
    plt.title('strain_lines, (eiv 3)')

    plt.subplot(322)
    plt.imshow(eigval1)
    for i in xrange(stretch_lines.shape[0]):
        plt.plot(stretch_lines[i, 1, :], stretch_lines[i, 0, :], 'w.', ms=1)
    plt.colorbar()
    plt.title('strain_lines, (eiv 1)')



    plt.subplot(323)
    plt.imshow(vel[:, :, 14, 0, 0])
    plt.colorbar()

    plt.subplot(324)
    plt.imshow(vel[:, :, 14, 1, 0])
    plt.colorbar()

    plt.subplot(325)
    plt.imshow((initpts3))
    plt.colorbar()
    plt.title('initpts3')

    plt.subplot(326)
    plt.imshow(ftle, vmin=0, vmax=8, cmap='jet', aspect='auto')
    plt.colorbar()
    plt.title('ftle')





    # plt.subplot(235)
    # plt.imshow(magU)
    # for j in xrange(strain_lines.shape[0]):
    #     plt.plot(ellipticp_lines[17, 1, :], strain_lines[17, 0, :], 'ro-', ms=1)
    #
    # plt.subplot(236)
    # plt.imshow(magU)
    # for j in xrange(strain_lines.shape[0]):
    #     plt.plot(ellipticm_lines[17, 1, :], strain_lines[17, 0, :], 'ro-', ms=1)

    plt.show()
