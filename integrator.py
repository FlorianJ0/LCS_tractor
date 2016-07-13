from scipy import ndimage

import numpy as np
from matplotlib import pyplot as plt


# sys.path.append('pytricubic-master/')
# import tricubic

def Integrator(vel, z, tphys, dt, nx, ny, nz, dim, domain, simtstep):
    print(z, tphys, dt, nx, dim, domain, simtstep)
    print(vel.shape)
    coeff = 1  # multiplicateur du nb de points
    # xi, yi, zi = np.mgrid[domain[0]:domain[1]:50*1j, domain[2]:domain[3]:50*1j, domain[4]:domain[5]:50*1j]
    xi, yi, zi = np.mgrid[0:nx:nx * coeff * 1j, 0:ny:ny * coeff * 1j, 0:nz:nz * coeff * 1j]
    positions = np.vstack([xi.ravel(), yi.ravel(), zi.ravel()]).T
    print positions.shape

    def vect(coord, t):
        x = [coord[0]]
        y = [coord[1]]
        z = [coord[2]]
        ux = ndimage.map_coordinates(vel[:, :, :, 0, :], [x, y, z, [t]], order=3)
        uy = ndimage.map_coordinates(vel[:, :, :, 1, :], [x, y, z, [t]], order=3)
        uz = ndimage.map_coordinates(vel[:, :, :, 2, :], [x, y, z, [t]], order=3)
        return [ux[0], uy[0], uz[0]]

    print vect([25, 37, 21], 2)
    print vel[25, 37, 22, :, 2]
    dt = 0.1
    t0 = 0
    t1 = vel.shape[-1]
    t = np.arange(t0, t1 + dt, dt)
    print t.shape

    start_pts = positions[62500:62525]
    print start_pts.shape
    trajectories = []
    trajectories1 = []

    def euler(f, x0, t0, t1, n):
        h = (t1-t0)/n
        t = np.arange(t0, t1+h, h)
        # print t,n , len(t)
        # return
        x = np.zeros([n,3])
        for i in xrange(int(n)-1):
            t[i+1] = t[i] + h
            x[i+1] = np.array(f(x[i],t[i])) * h

        return (t,x)


    for start_pt in start_pts:
        print 'prout'
        # trajectory = odeint(vect, start_pt, t)
        trajectory = euler(vect, start_pt, 2, vel.shape[-1], 20.) #FLOAT !!!
        trajectories.append(trajectory)

    fig = plt.figure()
    # ax = fig.add_subplot(1, 2, 1, projection='3d')
    # ax.plot(x, y, z)

    ax = fig.add_subplot(1, 1, 1, projection='3d')
    for trajectory in trajectories:
        x, y, z = trajectory.T  # transpose and unpack
        # x, y, z = zip(*trajectory)  # this also works!
        ax.plot(x, y, z)

    plt.show()
