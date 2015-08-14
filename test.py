# import numpy as np
# tt=5
# dim=2
# xmax = 3
# ymax = 3
# # int step
# h = 0.1
# # generate flow field
# # vx = np.empty((3,3,tt))
# # vy = np.empty((3,3,tt))
# # vz = np.empty((3,3,tt))
#
# x = np.indices((3,3))
#
# U = np.empty((dim,xmax,ymax,tt))
# # vectlist = [vx, vy]
# # print vectlist
# for t in xrange(tt):
#     for i in xrange(2):
#         toto=np.random.rand(9)
#         U[i,:,:,t]=toto.reshape((3,3))
#
#
# # initial point list
# # ptlist = x
# newptlist = x
# X=0
# # print newptlist[0]
# # # U=[vx, vy, vz]
# # # print U[0]
# # print newptlist.shape
#
# # print newptlist[0,:,:]
# # print U[0,:,:,0]
# for t in xrange(tt):
#     print 'time=', t
#     for i in xrange(dim):
#         newptlist[i,:,:] =np.multiply(newptlist[i,:,:],U[i,:,:,t] ) *h
#
# # print X
# import pygsl.sf
# # print "%g+/-%g"%pygsl.sf.erf(1)
import numpy as np

from vtk import vtkRectilinearGridReader
from vtk.util import numpy_support as vn

reader = vtkRectilinearGridReader()
reader.SetFileName('data/visit_ex_db00001.vtk')
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
# print 'file:', i
print 'dim:', dim
arrowglyph = data.GetPointData().GetArray('internalMesh/U')
vectU = vn.vtk_to_numpy(arrowglyph)
vectU[vectU < -20] = 'Nan'
U = np.empty((nx, ny, nz, 3, 1))
count = 0
s=0
i=j=k=0
for i in xrange(nz):
    for j in xrange(ny):
        for k in xrange(nz):
            U[i, j, k, :, count] = vectU[s]
            s+=1
count += 1