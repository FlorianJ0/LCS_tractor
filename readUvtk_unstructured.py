__author__ = 'p0054421'

# import os
import numpy as np

# from vtk import vtkRectilinearGridReader
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as vn
from scipy.interpolate import LinearNDInterpolator
import psutil
# import
# import shutil
import glob
import tricubic



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

        def swaap():
            if psutil.virtual_memory()[2]>80:
                print 'swaaaaaaaaaaaaap'
                quit()

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

        grid_x, grid_y, grid_z = np.mgrid[xmin:xmax:nx*1j,ymin:ymax:ny*1j,zmin:zmax:nz*1j]
        toto = np.array([pos,vectU])
        print "b4", toto.nbytes
        swaap()
        k=0
        for k in xrange(vectU.shape[0]):
            if toto[0,k,2] < zmin or toto[0,k,2] > zmax:
                toto[:,k,:] = np.nan
            k+=1
        # a[~np.isnan(a).any(axis=1)]
        toto=toto[~np.isnan(toto).any(axis=2)]
        print 'plein de nan ajoutes'
        print "after", toto.nbytes

        swaap()
        # grid_z1 = interpn(toto[0,:,:], toto[1,:,0], (grid_x, grid_y, grid_z), method='linear')
        # ex = LinearNDInterpolator((x, y, z), v)


        quit()
        # U=np.flipud(U)
        count += 1
    dt = 0.04
    tphys = count * dt
    print 'tphys = ', tphys
    print 'Files read'
    print '<<<<<<<<<>>>>>>>>>>'

    return