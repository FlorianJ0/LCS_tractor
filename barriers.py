__author__ = 'p0054421'
import numpy as np
import time
import matplotlib.pyplot as plt


def helic(vect, dx, dy, dz):
    print 'vect shape', vect.shape
    nx = vect.shape[0]
    ny = vect.shape[1]
    dd=np.empty([nx,ny,3])
    ddd=np.empty([nx,ny])
    eps0 = 1e-2
    # at each point, 3d grad vector of 3D eigenc
    for i in xrange(nx):
        for j in xrange(ny):
            dd[i,j,:]=np.gradient(vect[i,j,:])
            ddd[i,j]=np.dot(vect[i,j,:],dd[i,j,:])
            # if abs(ddd[i,j])<eps0:
            #     ddd[i,j]=True
            # else
            #     ddd[i,j]=False
    a = abs(ddd)<eps0
    ddd *= a

    return ddd

def reduced_lines(vect, dx, dy, dz, initpts):
    print 'integrate reduced LCSs'
    #on suppose qu on est toujours normal  a z
    # norm vect = 0 0 -1
    # donc n vectproduct k = kj -ki 0
    # la trajectoire est portee par le vect kj -ki 0 donc dans le plan
    # def func:

    return


def barrier_type(toto, eigval1, eigval3, eigvec1, eigvec3, tphys, dt, nx, ny, nz, domain, simtstep):
    dx = abs(domain[1] - domain[0]) / nx
    dy = abs(domain[3] - domain[2]) / ny
    dz = abs(domain[5] - domain[4]) / ny
    if toto == 0:
        print 'hyperbolic repelling LCS'
        # cond helicity
        # H = np.dot(eigvec[0],np.cross(np.gradient() ))
        # eigvecm = eigvec[:,:,0,:]
        # print eigvec.shape
        # print eigvecm.shape
        initpts = helic(eigvec3, dx, dy, dz)
        # f, ((ax1, ax2)) = plt.subplots(2)
        uu = initpts
        print uu.nonzero()
        print '000'
        print uu[uu.nonzero()]
        print '000'
        print uu.nonzero().shape, uu[uu.nonzero()].shape
        # imgplot = plt.imshow(halp)
        f, (ax1, ax2) = plt.subplots(2, 1,sharey=True)
        ax1.imshow(np.abs(uu))
        ax2.imshow(uu)
        plt.show()
        reduced_lines(eigvec3, dx, dy, dz, initpts)
        # r = eigenVectors[0]
    elif toto == 1:
        print 'hyperbolic attracting LCS'
        halp = helic(eigvec1, dx, dy, dz)
        # f, ((ax1, ax2)) = plt.subplots(2)
        uu = halp
        imgplot = plt.imshow(halp)
    elif toto == 2:
        print 'ellipc LCS'
    else:
        print 'tu sais pas coder'
        print 'error', np.random.random_integers(0, 100)
