import numpy as np
# from vtk.util import numpy_support as vn
import time
from vtk.util.numpy_support import vtk_to_numpy
# import shutil
from vtk import *
import vtk
from tempfile import TemporaryFile

outfile = TemporaryFile()
fname = '/home/p0054421/MEGA/calcul/LCS_tractor/data/paths.vtk'
npArray = '/home/p0054421/MEGA/calcul/LCS_tractor/data/paths'
# fname = '/home/florian/MEGAsync/calcul/LCS_tractor/data/paths.vtk'
if os.path.isfile(npArray + '.npy'):
    print 'loading ', npArray
    t0 = time.time()
    a = np.load(npArray + '.npy')
    print 'already existing array loaded in', time.time() - t0, 's'
else:
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fname)
    reader.Update()
    data = reader.GetOutput()
    ll = data.GetLines()
    n_pts = ll.GetSize()  # nb of points
    n_lines = ll.GetNumberOfCells()  # nb of lines
    idList = vtk.vtkIdList()
    idList.SetNumberOfIds(n_lines)
    print idList
    idList.SetId(0, 0)
    abscissaArray = vtk.vtkFloatArray()
    traj = []  # list of pos
    vel = []  # list of vel
    vort = []  # list of vort
    k = 0
    for i in xrange(n_lines):
        cell = data.GetCell(i)
        abscissa = 0.0
        previousPoint = None
        if i % np.int(n_lines / 100) == 0:
            print k, '% read'
            k += 1
        print 'shit\'s read'

        llength = cell.GetNumberOfPoints()
        pos = []  # np.empty([llength, 3])  # pos, u, vort
        u = []  # np.empty([llength, 3])  # pos, u, vort
        vorti = []  # np.empty([llength])  # pos, u, vort
        for j in range(llength):
            pointId = cell.GetPointId(j)
            pos.append(data.GetPoint(pointId))
            u.append(vtk_to_numpy(data.GetPointData().GetArray("U"))[pointId])
            vorti.append(vtk_to_numpy(data.GetPointData().GetArray("Vorticity"))[pointId])

        traj.append(pos)
        vel.append(u)
        vort.append(vorti)

    x = [[traj], [vel], [pos]]
    xx = np.array(x)
    np.save(npArray, x)


# print idList.SetId(3)
# print n_pts, n_lines
