__author__ = 'p0054421'
import numpy as np
import time
import matplotlib.pyplot as plt
from airkaeffe import rk45, heun, euler



def helic(vect, dx, dy, dz):
    print 'vect shape', vect.shape
    nx = vect.shape[0]
    ny = vect.shape[1]
    dd = np.empty([nx, ny, 3])
    ddd = np.empty([nx, ny])
    eps0 = 1e-2
    # at each point, 3d grad vector of 3D eigenc
    for i in xrange(nx):
        for j in xrange(ny):
            dd[i, j, :] = np.gradient(vect[i, j, :])
            ddd[i, j] = np.dot(vect[i, j, :], dd[i, j, :])
            # if abs(ddd[i,j])<eps0:
            #     ddd[i,j]=True
            # else
            #     ddd[i,j]=False
    a = abs(ddd) < eps0
    ddd *= a

    return ddd



def reduced_lines(vect, nx, ny, nz, initpts):
    def gamma(x, t):
        newpos = np.empty((2))
        cond = np.dot(vect[x[0], x[1], :], np.gradient(vect[x[0], x[1], :]))
        if (abs(cond) < 5e-1):
            #i did the maths, must be right
            newpos[0] = - 1. * vect[x[0], x[1], 1]
            newpos[1] =  1. * vect[x[0], x[1], 0]
            # print cond
            # print newpos
            # print - 1. * vect[x[0], x[1], 1]
            # print '000000000000000000000000000'

        else:
            newpos[0] = x[0]
            newpos[1] = x[1]
        # print 'cond', cond
        return newpos
    print 'integrate reduced LCSs'
    N = 50
    # on suppose qu on est toujours normal  a z
    # norm vect = 0 0 -1
    # donc n vectproduct k = kj -ki 0
    # la trajectoire est portee par le vect kj -ki 0 donc dans le plan
    t = np.linspace(0, 1, N) #pseudo integrator
    line = np.zeros ((initpts.shape[1], 2, N))

    for i in xrange(initpts.shape[1]):
        y0 = initpts[:, i]*1.
        line[i, :, :] = euler(gamma, y0, t).swapaxes(1,0)
        print y0
        print line[i, :, :]
        print '-*---------------*-'
        #
        # print i
        # print y0
        # print line
        # print line.shape
        # print line[-1]

    return line

def barrier_type(toto, eigval1, eigval3, eigvec1, eigvec3, tphys, dt, nx, ny, nz, domain, simtstep):
    dx = abs(domain[1] - domain[0]) / nx
    dy = abs(domain[3] - domain[2]) / ny
    dz = abs(domain[5] - domain[4]) / ny
    print 'hyperbolic  LCS'
    # cond helicity
    # H = np.dot(eigvec[0],np.cross(np.gradient() ))
    # eigvecm = eigvec[:,:,0,:]
    # print eigvec.shape
    # print eigvecm.shape
    initpts = helic(eigvec3, dx, dy, dz)
    seeds = np.array([initpts.nonzero()[0], initpts.nonzero()[1]])
    strain_lines = reduced_lines(eigvec3, nx, ny, nz, seeds)
    initpts0 = helic(eigvec1, dx, dy, dz)
    seeds0 = np.array([initpts0.nonzero()[0], initpts0.nonzero()[1]])
    stretch_lines = reduced_lines(eigvec1, nx, ny, nz, seeds0)
    print stretch_lines.shape
    k = 0
    toto=np.empty((stretch_lines.shape[0]*stretch_lines.shape[1],2))
    for i in xrange(stretch_lines.shape[0]):
        for j in xrange(stretch_lines.shape[1]):
            toto[k,:] = stretch_lines[i,:,k]
            k+=1

    f, (ax1, ax2) = plt.subplots(2, 1, sharey=True)
    ax1.imshow(np.abs(initpts))
    ax2.imshow(toto)
    plt.show()


