__author__ = 'p0054421'
import numpy as np
import time
import matplotlib.pyplot as plt


def helic(vect, dx, dy, dz):
    print 'vect shape', vect.shape
    nx = vect.shape[0]
    ny = vect.shape[1]
    dd=np.empty([nx,ny,3])
    ddd=np.empty([nx,ny,3])
    # at each point, 3d grad vector of 3D eigenc
    for i in xrange(nx):
        for j in xrange(ny):
            dd[i,j,:]=np.gradient(vect[i,j,:])
            ddd[i,j]=np.dot(vect[i,j,:],dd[i,j,:])
    return ddd



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
        halp = helic(eigvec3, dx, dy, dz)
        # f, ((ax1, ax2)) = plt.subplots(2)
        uu = halp
        imgplot = plt.imshow(halp)
        # r = eigenVectors[0]
    elif toto == 1:
        print 'hyperbolic attracting LCS'
    elif toto == 2:
        print 'ellipc LCS'
    else:
        print 'tu sais pas coder'
        print 'error', np.random.random_integers(0, 100)
