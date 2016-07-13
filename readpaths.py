__author__ = 'p0054421'

# import os
import numpy as np
from vtk import vtkPolyDataReader
# from vtk.util import numpy_support as vn
from vtk.util.numpy_support import vtk_to_numpy
# import shutil
import glob


def read_files(loc):
    """

    :rtype : vector
    """
    dircont = sorted(glob.glob(loc + '/*.vtk'))
    pat = sorted(glob.glob(loc + '/paths.vtk'))
    if len(pat) == 0:
        print 'po de fichiers'
        print 'con'
        #quit()
    count = 0

    for i in pat:
        reader = vtkPolyDataReader()
        # reader = vtkStructuredPointsReader()
        reader.SetFileName(i)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        uu =reader.GetOutput().GetPointData().GetArray("U")
        U = vtk_to_numpy(uu)
        vv = reader.GetOutput().GetPointData().GetArray("Vorticity")
        # vtk_to_numpy(data.GetPointData().GetArray('p'))[1000]

        vort = vtk_to_numpy(vv)
        print U.shape
        print vort.shape

        quit()
        # dim = data.GetDimensions()

        # nx = dim[0] #- 1
        # ny = dim[1] #- 1
        # nz = dim[2] #- 1
        #
        # bounds = data.GetBounds()
        # xmin = bounds[0]
        # xmax = bounds[1]
        # ymin = bounds[2]
        # ymax = bounds[3]
        # zmin = bounds[4]
        # zmax = bounds[5]
        # domain = np.array([xmin, xmax, ymin, ymax, zmin, zmax])

        print 'file:', i
        # print 'dim:', dim

        # arrowglyph = data.GetCellData().GetArray('internalMesh/U')
        arrowglyph = data.GetPointData().GetArray('U')

        # print data.GetCellData()
        vectU = vn.vtk_to_numpy(arrowglyph)

        # vectU[vectU < -20] = 'NaN'
        # if np.all(vectU) < np.array([1e-8,1e-8,1e-8])):
        if (count == 0):
            U = np.empty((nx, ny, nz, len(dim), len(dircont)))

        s = i = j = k = 0

        for k in xrange(nz):
            for j in xrange(ny):
                for i in xrange(nx):
                    U[i, j, k, :, count] = vectU[s]
                    s += 1
        count += 1

    dt = 0.05
    tphys = count * dt
    print 'tphys = ', tphys
    print 'Files read'
    print '<<<<<<<<<>>>>>>>>>>'

    return domain
