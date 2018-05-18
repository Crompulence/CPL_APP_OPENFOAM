import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os

#Run simulation
fDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(fDir)
os.system('./run.sh')

#Add pyDataView to system path
ppdir = '/home/asufian/Programs/coupled_LAMMPS_OpenFOAM/pyDataView'
sys.path.append(ppdir)
import postproclib as ppl

#Get post-processing object
OPEN_FOAM_CASE = './openfoam/'
PPObj = ppl.OpenFOAM_PostProc(OPEN_FOAM_CASE)

#Extract the pressure field object
pObj = PPObj.plotlist['p']

#Extract grid configuration and physical size in y-direction (gravity direction)
nx = pObj.Raw.ncx
ny = pObj.Raw.ncy
nz = pObj.Raw.ncz
ylo = 0
yhi = pObj.Raw.yL

#Extract pressure field for the final timestep
p = pObj.read(startrec = pObj.maxrec, endrec = pObj.maxrec)
p = np.reshape(p, [nx, ny, nz], order='F')
pMean = np.mean(np.mean(p, 2), 0)

#Solution to hydrostatic pressure field
with open(OPEN_FOAM_CASE + 'constant/environmentalProperties', 'r') as fObj:
    for l in fObj:
        if l.startswith('g'):
            g = l[l.find('(')+1:l.find(')')]
            g = [float(i) for i in g.split()]
            gIdx = [i for i, x in enumerate(g) if (x != 0)][0]
            g = g[gIdx]

with open(OPEN_FOAM_CASE + 'constant/transportProperties', 'r') as fObj:
    for l in fObj:
        if 'rhob' in l:
            rhof = l.split()[-1]
            rhof = float(rhof[:-1])

h = np.linspace(yhi, ylo, ny+1)
h = 0.5*(h[1:ny+1] + h[0:ny])
pSol = np.tile(-rhof*g*h, [nx, nz, 1]);
pSol = np.swapaxes(pSol, 1, 2);
pSolMean = np.mean(np.mean(pSol, 2), 0)

#Check that numerical results match hydrostatic solution to some tolerance
def test_hydrostatic_complete():
	tol = 0.01
	err = abs((p - pSol)/pSol) <= tol
	assert(err.all())
	print(str(nx*ny*nz - np.count_nonzero(err)) + ' of ' + str(nx*ny*nz) + ' grid cells exceed the specified error tolerance of ' + str(tol*100) + '%')

def test_hydrostatic_average():
	tol = 0.01
	err = abs((pMean - pSolMean)/pSolMean) <= tol
	assert(err.all())
	print(str(ny - np.count_nonzero(err)) + ' of ' + str(ny) + ' averaged layers exceed the specified error tolerance of ' + str(tol*100) + '%')

#Plot comparison of numerical and analytical solution
plt.plot(h, pMean, 'r-', linewidth=3.0)
plt.plot(h, pSolMean, 'k--')
plt.xlabel('Depth (cm)')
plt.ylabel('Pressure (gcm$^{-1}$s$^{-2}$)')
plt.legend(('Numerical', 'Analytical'))
plt.savefig('fig_hydrostatic.png')
plt.close()