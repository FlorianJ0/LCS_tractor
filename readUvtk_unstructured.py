__author__ = 'p0054421'

# import os
import numpy as np

# from vtk import vtkRectilinearGridReader
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as vn
from scipy.interpolate import RegularGridInterpolator

# import shutil
import glob


def read_files(loc, dim, extend):
    """

    :rtype : vector
    """
    # file = '/media/backup/patients_article0/patient4/DOIRE^JEAN-LOUIS/DOIRE^JEAN_LOUIS_20060124/Simulation/VTK/Simulation_9322.vtk'
    # print file
    dircont = sorted(glob.glob(loc + '/*.vtk'))
    if len(dircont) == 0:
        print 'po de fichiers'
        print 'con'
        quit()
    count = 0
    for i in dircont:
        reader=vtkUnstructuredGridReader()
        reader.SetFileName(i)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        data = reader.GetOutput()
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
        arrowglyph = data.GetPointData().GetArray('U')
        coord = data.GetPoint
        vectU = vn.vtk_to_numpy(arrowglyph)

        pos = np.empty((vectU.shape[0], 3))
        for l in xrange(vectU.shape[0]):
            pos[l][:] = np.array(coord(l))

        if pos.shape != vectU.shape:
            print 'oups !'
            quit()

        print count
        # U=np.flipud(U)
        count += 1
    dt = 0.04
    tphys = count * dt
    print 'tphys = ', tphys
    print 'Files read'
    print '<<<<<<<<<>>>>>>>>>>'

    return