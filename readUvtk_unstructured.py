__author__ = 'p0054421'

# import os
import numpy as np

# from vtk import vtkRectilinearGridReader
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as vn
from scipy.interpolate import griddata
import psutil
import
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
        # bounds = data.GetBounds()
        xmin = extend[0]
        xmax = extend[1]
        ymin = extend[2]
        ymax = extend[3]
        zmin = extend[4]
        zmax = extend[5]
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
        if psutil.swap_memory()[3]>10:
            print 'swaaaaaaaaaaaaap'
            quit()
        grid_x, grid_y, grid_z = np.mgrid[xmin:xmax:nx*1j,ymin:ymax:ny*1j,zmin:zmax:nz*1j]

        grid_z1 = griddata(pos, vectU, (grid_x, grid_y, grid_z), method='linear')
        # U=np.flipud(U)
        count += 1
    dt = 0.04
    tphys = count * dt
    print 'tphys = ', tphys
    print 'Files read'
    print '<<<<<<<<<>>>>>>>>>>'

    return