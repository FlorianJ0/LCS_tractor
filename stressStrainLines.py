import numpy as np
import sympy.physics.vector as sp



def reducedLineCriteria(vectSlice, eps=1e-4):
    #computes the starting points for the strain/stress Lines
    #pts must satisfy abs(H(ki3) = abs(avg(grad(ki3) cross ki3)) < eps
    #H  = helicity
    Helic = np.dot(vectSlice,sp.curl(vectSlice))
    return


def strainLines(vel, CGtensor, eigval1, eigval3, eigvec1, eigvec3, domain, tphys, dt,simtstep):
    print "going for strain lines"
    bobol = vel[np.abs(vel[:, :, :, 2, 0])<1e-6]
    nx, ny, nz = vel.shape[:3]
    print nx, ny, nz
    reducedLineCriteria(eigvec3[:,:,12,:], eps=1e-4)
    return