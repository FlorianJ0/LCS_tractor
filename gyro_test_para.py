__author__ = 'p0054421'

__author__ = 'p0054421'

import numpy as np
import time
import sys
import pp


nx=50
ny =25
nz = 5
dx = 4e-2
dy = 4e-2
dz = 4e-2
ttot = 200
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
    stamp=time.time()
    vel = np.empty((nx,ny,2,int(ttot/dt))) #
    nz = 5
    vel3D = np.empty((nx,ny, nz, 3,int(ttot/dt))) #

    for t in range(0, int(ttot/dt)):
        tt=t*dt
        # print tt
        for i in xrange(nx):
            for j in xrange(ny):
                x = i * dx
                y = j * dy

                vel[i,j,0,t] = - np.pi*0.1*np.sin(np.pi*f(x, tt))*np.cos(np.pi*y)*1
                vel[i,j,1,t] =   np.pi*0.1*np.cos(np.pi*f(x, tt))*np.sin(np.pi*y)*dfdx(x, tt)*1
        vel3D[:,:,0,0:2,t] = vel[:,:,:,t]

    for k in range(1, nz):
        vel3D[:,:,k,0:2,t] = vel3D[:,:,0,0:2,t]
    vel3D[:,:,:,2,:] = -1e-3
    time_tot=time.time()-stamp
    print 'time_tot seq %f s ' %time_tot
    return time_tot

def gyro_par():
    nx=50
    ny =25
    nz = 5
    dx = 4e-2
    dy = 4e-2
    dz = 4e-2
    ttot = 200
    dt = 0.5
    w=2*np.pi/10
    stamp=time.time()
    vel = np.empty((nx,ny,2,int(ttot/dt))) #
    nz = 5
    vel3D = np.empty((nx,ny, nz, 3,int(ttot/dt))) #

    for t in range(0, int(ttot/dt)):
        tt=t*dt
        # print tt
        for i in xrange(nx):
            for j in xrange(ny):
                x = i * dx
                y = j * dy

                vel[i,j,0,t] = - np.pi*0.1*np.sin(np.pi*f(x, tt))*np.cos(np.pi*y)*1
                vel[i,j,1,t] =   np.pi*0.1*np.cos(np.pi*f(x, tt))*np.sin(np.pi*y)*dfdx(x, tt)*1
        vel3D[:,:,0,0:2,t] = vel[:,:,:,t]

    for k in range(1, nz):
        vel3D[:,:,k,0:2,t] = vel3D[:,:,0,0:2,t]
    vel3D[:,:,:,2,:] = -1e-3
    time_tot=time.time()-stamp
    print 'time_tot para%f s ' %time_tot

    return time_tot


toto = gyro()
# print 'time sep %f s' %toto

ppservers = ()

job_server = pp.Server(8, ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"
# gyro()
job_server.submit(gyro_par(),args=(), depfuncs=(), modules=("numpy as np", "time"),callback=None, callbackargs=(), group='default', globals=None)
# toto = job1()
# print job1
