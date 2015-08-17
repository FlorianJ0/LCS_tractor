__author__ = 'p0054421'

__date__ = '30 Juillet 2015'

'obj: calculer les  LCS hyperboliques depuis un champ de vitesse, ecoulement incompressible, periodique de sang dans un' \
'AAA. Si possible ajouter LCS  elliptic et parabolic si j ai le temps. Le tout sur du vtk structured'
#
import readUvtk
import cauchygreen
import barriers
import cauchygreen3d
import gyronator
# vitesse min pour savoir si on est in/out domain (a mettre a -100 pour la suite)
outofdomain = 1E-16
instat = True
# t step de la simu
simtstep = 0.04
# read file/s
loc = '/home/p0054421/MEGA/calcul/LCS_tractor/data'
vel, nx, ny, nz, dim_initial, tphys, dt, domain = readUvtk.read_files(loc)
# vel, nx, ny, nz, dim_initial, tphys, dt, domain = gyronator.gyro()
print 'Velocity read'

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
z = 25.
# t = 3
# calcul sur un plan x y parceque jsuis trop une feignasse pour un code generique
dim = 2
# eigval1, eigval3, eigvec1, eigvec3, interpU_i = cauchygreen3d.cgstki3(vel, z, tphys, dt, nx, ny, nz, 3, domain, simtstep)
eigval1, eigval3, eigvec1, eigvec3, interpU_i = cauchygreen3d.cgstki3(vel, z, tphys, dt, nx, ny, nz, 3, domain, simtstep)

barriers.barrier_type(0, eigval1, eigval3, eigvec1, eigvec3, tphys, dt, nx, ny, nz, domain, simtstep)
# barriers.barrier_type(0, eigval, eigvec, tphys, dt, nx, ny, domain, simtstep)
