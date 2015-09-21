__author__ = 'p0054421'
import time
import ConfigParser

import numpy as np
import matplotlib.pyplot as plt



# from Scientific.Functions.Interpolation import InterpolatingFunction as IF
from scipy import ndimage

from airkaeffe import heun

Config = ConfigParser.ConfigParser()
Config.read('parameters.ini')


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print "skip: %s" % option
        except:
            print "exception on %s!" % option
            dict1[option] = None
    return dict1


def rpog():
    print("# "),


def helic(vect, dx, dy, dz, eps0, bobol):
    print 'vect shape', vect.shape
    nx = vect.shape[0]
    ny = vect.shape[1]
    dd = np.empty([nx, ny, 3])
    ddd = np.empty([nx, ny])
    # eps0 = 5e-1
    # at each point, 3d grad vector of 3D eigenc
    for i in xrange(nx):
        for j in xrange(ny):
            if bobol[i, j]:
                # dd[i, j, :] = np.gradient(vect[i, j, :])
                ddd[i, j] = np.dot(vect[i, j, :], np.gradient(vect[i, j, :]))
            else:
                ddd[i, j] = 0
                # if abs(ddd[i,j])<eps0:
                #     ddd[i,j]=True
                # else
                #     ddd[i,j]=False

    npp = ConfigSectionMap('cauchygreen')['nseeds']
    # print'sort'
    # print ddd.shape
    # sort = insertion_sort(ddd)[npp]
    # print sort
    a = abs(ddd) < eps0
    ddd *= a

    return ddd


def reduced_lines(vect, nx, ny, nz, initpts, thresh, n, bobol):
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)
    vecti = (vect[:, :, 0])
    # vecti /=np.max(np.abs(vecti))
    vectj = (vect[:, :, 1])
    # vectj /=np.max(np.abs(vectj))

    # vectj = vect[:,:,1]
    axes = (x, y)
    # vi = IF(axes, vecti)
    # vj = IF(axes, vectj)
    # vi = ndimage.map_coordinates(vecti)

    # vj = IF(axes, vectj)


    def gamma(x, t):
        newpos = np.empty(2)
        # print 'x', x
        vi = ndimage.map_coordinates(vecti, [[x[0]], [x[1]]], order=2, mode='constant', cval=0.0, prefilter=False)
        vj = ndimage.map_coordinates(vectj, [[x[0]], [x[1]]], order=2, mode='constant', cval=0.0, prefilter=False)

        # a = np.array([vi(x[0], x[1]), vj(x[0], x[1])])
        a = np.array([vi, vj])[:, 0]

        cond = np.dot(a, np.gradient(a))
        cut = False
        # if (1e-4 < abs(cond) < thresh):
        if abs(cond) < thresh and bobol[x[0], x[1]]:
            # i did the maths, must be right
            newpos[0] = - 1. * vj  # (x[0], x[1])
            newpos[1] = 1. * vi  # (x[0], x[1])
        else:
            cut = True
        # print newpos
        return newpos, cut

    print 'integrate reduced LCSs'
    N = n * 1
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
        # print 'line number %i' % i
        # print y0
        # print line[i, :, -1]
        # print np.sqrt((np.abs(y0[0]-line[i, 0, -1]))**2+(np.abs(y0[1]-line[i, 1, -1]))**2)
        # print line.shape
        # fname = 'toto' + str(i)
        # np.savetxt(fname, line[i,:,:])
        a = 1. * i / initpts.shape[1]
        if (100 * a) % 10 < 1e-1:
            rpog()

    print '\n'
    return line


def barrier_type(toto, eigval1, eigval3, eigvec1, eigvec3, vel, tphys, dt, nx, ny, nz, domain, simtstep, bobol):
    dx = abs(domain[1] - domain[0]) / nx
    dy = abs(domain[3] - domain[2]) / ny
    dz = abs(domain[5] - domain[4]) / ny
    # print 'hyperbolic  LCS'
    # cond helicity
    # H = np.dot(eigvec[0],np.cross(np.gradient() ))
    # eigvecm = eigvec[:,:,0,:]
    # print eigvec.shape
    # print eigvecm.shape

    a = np.sqrt(np.sqrt(eigval1)) / (np.sqrt(eigval1) + np.sqrt(eigval3))
    b = np.sqrt(np.sqrt(eigval3)) / (np.sqrt(eigval1) + np.sqrt(eigval3))
    nnp = np.empty((nx, ny, 3))
    nnm = np.empty((nx, ny, 3))
    # print a.shape, nnp.shape
    for i in range(3):
        nnp[:, :, i] = a * eigvec1[:, :, i] + b * eigvec3[:, :, i]
        nnm[:, :, i] = a * eigvec1[:, :, i] - b * eigvec3[:, :, i]

    # threshhods of seed finding
    ths_strain_lines = np.float(ConfigSectionMap('barriers')['ths_strain_lines'])
    ths_stretch_lines = np.float(ConfigSectionMap('barriers')['ths_stretch_lines'])
    ths_ellipticp_lines = np.float(ConfigSectionMap('barriers')['ths_ellipticp_lines'])
    ths_ellipticm_lines = np.float(ConfigSectionMap('barriers')['ths_ellipticn_lines'])

    # threshhods of line integration
    th_strain_lines = np.float(ConfigSectionMap('barriers')['th_strain_lines'])
    th_stretch_lines = np.float(ConfigSectionMap('barriers')['th_stretch_lines'])
    th_ellipticp_lines = np.float(ConfigSectionMap('barriers')['th_ellipticp_lines'])
    th_ellipticn_lines = np.float(ConfigSectionMap('barriers')['th_ellipticn_lines'])

    stamp = time.time()
    initpts3 = helic(eigvec3, dx, dy, dz, ths_strain_lines, bobol)
    # print initpts3.nonzero()
    seeds3 = np.array([initpts3.nonzero()[0], initpts3.nonzero()[1]])
    print seeds3.shape
    strain_lines = reduced_lines(eigvec3, nx, ny, nz, seeds3, th_strain_lines, 500, bobol)
    print "number of seeds %i" % seeds3.shape[1]
    print '-----------------------------------------------------'
    print 'strain lines (repelling) computed  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'

    stamp = time.time()
    initpts1 = helic(eigvec1, dx, dy, dz, ths_stretch_lines, bobol)
    seeds1 = np.array([initpts1.nonzero()[0], initpts1.nonzero()[1]])
    stretch_lines = reduced_lines(eigvec1, nx, ny, nz, seeds1, th_stretch_lines, 500, bobol)
    print 'number of seeds %i' % seeds1.shape[1]
    print '-----------------------------------------------------'
    print 'stretch (attracting) lines computed  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'

    stamp = time.time()
    initptsp = helic(nnp, dx, dy, dz, ths_ellipticp_lines, bobol)
    seedsp = np.array([initptsp.nonzero()[0], initptsp.nonzero()[1]])
    ellipticp = reduced_lines(nnp, nx, ny, nz, seedsp, th_ellipticp_lines, 500, bobol)
    print 'number of seeds %i' % seedsp.shape[1]
    print '-----------------------------------------------------'
    print 'ellipticp lines computed  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'

    stamp = time.time()
    initptsm = helic(nnm, dx, dy, dz, ths_ellipticm_lines, bobol)
    seedsm = np.array([initptsm.nonzero()[0], initptsm.nonzero()[1]])
    ellipticm = reduced_lines(nnm, nx, ny, nz, seedsm, th_ellipticn_lines, 500, bobol)
    print 'number of seeds %i' % seedsm.shape[1]
    print '-----------------------------------------------------'
    print 'ellipticp lines computed  in %f s ' % (time.time() - stamp)
    print '-----------------------------------------------------'

    ftle = (1 / 0.16) * np.log(np.sqrt(eigval3 + 1e-8))
    magU = np.sqrt(vel[:, :, 7, 0, 0] ** 2 + vel[:, :, 7, 1, 0] ** 2 + vel[:, :, 7, 2, 0] ** 2)

    plt.subplot(221)
    # eigval3 /=np.max(np.abs(eigval3-100))
    # eigval1 /=np.max(np.abs(eigval1))
    plt.imshow(magU)
    for j in xrange(strain_lines.shape[0]):
        plt.plot(strain_lines[:, 1, :], strain_lines[:, 0, :], 'k.', ms=1)
    plt.colorbar()
    plt.title('strain_lines, (eiv 3)')

    plt.subplot(222)
    plt.imshow(magU)
    for i in xrange(stretch_lines.shape[0]):
        plt.plot(stretch_lines[i, 1, :], stretch_lines[i, 0, :], 'k.', ms=1)
    plt.colorbar()
    plt.title('strain_lines, (eiv 1)')

    plt.subplot(223)
    plt.imshow(magU)
    for i in xrange(ellipticp.shape[0]):
        plt.plot(ellipticp[i, 1, :], ellipticp[i, 0, :], 'k.', ms=1)
    plt.colorbar()
    plt.title('ellipticp ')

    plt.subplot(224)
    plt.imshow(magU)
    for i in xrange(ellipticm.shape[0]):
        plt.plot(ellipticm[i, 1, :], ellipticm[i, 0, :], 'k.', ms=1)
    plt.colorbar()
    plt.title('ellipticm')

    #
    # plt.subplot(323)
    # plt.imshow(vel[:, :, 14, 0, 0])
    # plt.colorbar()
    #
    # plt.subplot(324)
    # plt.imshow(vel[:, :, 14, 1, 0])
    # plt.colorbar()
    #
    # plt.subplot(325)
    # plt.imshow((initpts3))
    # plt.colorbar()
    # plt.title('initpts3')
    #
    # plt.subplot(326)
    # plt.imshow(ftle, vmin=0, vmax=8, cmap='jet')
    # plt.colorbar()
    # plt.title('ftle')

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
    return
