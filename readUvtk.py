__author__ = 'p0054421'

# import os
import numpy as np

from vtk import vtkRectilinearGridReader
from vtk.util import numpy_support as vn

# import shutil
import glob


def read_files(loc):
    """

    :rtype : vector
    """
    dircont = sorted(glob.glob(loc + '/*.vtk'))
    count = 0
    for i in dircont:
        reader = vtkRectilinearGridReader()
        reader.SetFileName(i)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        data = reader.GetOutput()
        dim = data.GetDimensions()
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        bounds = data.GetBounds()
        xmin = bounds[0]
        xmax = bounds[1]
        ymin = bounds[2]
        ymax = bounds[3]
        zmin = bounds[4]
        zmax = bounds[5]
        domain = np.array([xmin, xmax, ymin, ymax, zmin, zmax])

        print 'file:', i
        print 'dim:', dim
        arrowglyph = data.GetPointData().GetArray('internalMesh/U')
        vectU = vn.vtk_to_numpy(arrowglyph)
        # vectU[vectU < -20] = 'NaN'
        if (count == 0):
            U = np.empty((nx, ny, nz, len(dim), len(dircont)))
        s = 0
        for k in xrange(nz):
            for j in xrange(ny):
                for i in xrange(nx):
                    U[j, i, k, :, count] = vectU[s]
                    s += 1
        # print U[25,25,25,0:2,count]
        print count
        # U=np.flipud(U)
        count += 1
    dt = 0.04
    tphys = count * dt
    print 'tphys = ', tphys
    print 'Files read'
    print '<<<<<<<<<>>>>>>>>>>'

    return U, nx, ny, nz, dim, tphys, dt, domain