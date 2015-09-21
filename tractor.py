__author__ = 'p0054421'

__date__ = '30 Juillet 2015'

'obj: calculer les  LCS hyperboliques depuis un champ de vitesse, ecoulement incompressible, periodique de sang dans un' \
'AAA. Si possible ajouter LCS  elliptic et parabolic si j ai le temps. Le tout sur du vtk structured'
#
import time
import ConfigParser

import numpy as np

import readUvtk
import barriers
import cauchygreen3d
import gyronator
import readUvtk_unstructured
import sup3D

ttot = time.time()
Config = ConfigParser.ConfigParser()
Config.read('parameters.ini')


# a = np.arange(10e7)
# b = np.arange(0, 20e7, 2)
#
# tttt = time.time()
# ne.ncores = 4
# c = ne.evaluate("2*a+3*b")
# print (time.time()-tttt)
# tttt = time.time()
# c = 2*a+3*b
# print (time.time()-tttt)
# quit()
# print c
def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print "skip: %s" % option
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


# vitesse min pour savoir si on est in/out domain (a mettre a -100 pour la suite)
outofdomain = ConfigSectionMap('filereading')['outofdomain']
# instat = True
# t step de la simu
simtstep = ConfigSectionMap('filereading')['timestep']
# read file/s
loc = ConfigSectionMap('filereading')['folder']

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Goord morning %s' % sup3D.hello()
print 'Welcome in my lair.'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
stamp = time.time()
struct = ConfigSectionMap('filereading')['datastruct']
if struct == 'unstruct':
    extend = np.array([0.05, 0.15, 0.04, 0.14, 0.2, 0.3])
    dim = np.array([20, 20, 20])
    vel, nx, ny, nz, dim_initial, tphys, dt, domain = readUvtk_unstructured.read_files(loc, dim, extend)
    # quit()
elif struct == 'struct':
    vel, nx, ny, nz, dim_initial, tphys, dt, domain = readUvtk.read_files(loc)
elif struct == 'anal':
    vel, nx, ny, nz, dim_initial, tphys, dt, domain = gyronator.gyro()
else:
    print 'wut ?'
    #quit()

print '-----------------------------------------------------'
print 'Velocity read in %f s ' % (time.time() - stamp)
print '-----------------------------------------------------'

"""
-sur [t0, t0+T] et n grilles G0 de PI 2D uniformes recti sur z (parceque ca m'arrange)
-compute cauchy-green strain tensor Ct0T, ki3 (dominant eiv) et ki1 (tant qu'a faire)
-selection d'une grille  reduce, ici dxr=2*dx et dyr=2*dy G1
a partir des points de G1, calcul de gammaS1 d'apres conditions (voir h2)
-integrer ces strainlines tant que |Hki3|<eps
- filter strainlines (Hausdorff distance criteria)
-repeat on s1
"""

# compute CG-strain tensor + eig/eiv on z plane
z = ConfigSectionMap('cauchygreen')['z']
z = float(z)
# t = 3
# calcul sur un plan x y parceque jsuis trop une feignasse pour un code generique
dim = 2

eigval1, eigval3, eigvec1, eigvec3, interpU_i, bobol = cauchygreen3d.cgstki3(vel, z, tphys, dt, nx, ny, nz, 3, domain,
                                                                             simtstep)
barriers.barrier_type(0, eigval1, eigval3, eigvec1, eigvec3, interpU_i, tphys, dt, nx, ny, nz, domain, simtstep, bobol)
print '-----------------------------------------------------'
print 'full total time  in %f s ' % (time.time() - ttot)
print '-----------------------------------------------------'
# barriers.barrier_type(0, eigval, eigvec, tphys, dt, nx, ny, domain, simtstep)
