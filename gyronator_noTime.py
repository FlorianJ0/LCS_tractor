__author__ = 'p0054421'

import numpy as np

nx=50
ny = 25
dx = 4e-2
dy = 4e-2
ttot = 15
dt = 0.5
w=2*np.pi/10
# w=0
def a(t):
    r =  0.5*np.sin(w*t)
    return r

def b(t):
    s = 1 - 2 * 0.5 * np.sin(w*t)
    return s

def f(x, t):
    toto = a(t) * x**2 + b(t)*x
    return toto

def dfdx(x, t):
    toto = 2 * a(t) * x + b(t)
    return toto



def gyro():
    vel = np.empty((nx,ny,2,ttot)) #
    for t in range(0, ttot):
        for i in xrange(nx):
            for j in xrange(ny):
                x = i * dx
                y = j * dy
                vel[i,j,0,t] = - np.pi*0.1*np.sin(np.pi*f(x, t))*np.cos(np.pi*y)
                vel[i,j,1,t] =   np.pi*0.1*np.cos(np.pi*f(x, t))*np.sin(np.pi*y)*dfdx(x, t)





    print vel.shape
    domain=np.array([0,2,0,1,1,1])
    nz = 0
    dim_initial = 2




    return vel, nx, ny, nz, dim_initial, ttot, domain
