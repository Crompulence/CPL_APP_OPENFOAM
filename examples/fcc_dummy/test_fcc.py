import numpy as np 
import sys
import os
import subprocess as sp

#Run simulation
fDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(fDir)
os.system('./run.sh')

#Add pyDataView to system path
sys.path.insert(0, "./pyDataView/")
try:
    import postproclib as ppl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/pyDataView.git ./pyDataView"
    downloadout = sp.check_output(cmd, shell=True)
    print(downloadout)
    sys.path.insert(0, "./pyDataView")
    import postproclib as ppl

#Get Post Proc Object
fdir = './openfoam/'
PPObj = ppl.OpenFOAM_PostProc(fdir)

#Get plotting object
pObj = PPObj.plotlist['p']
UbObj = PPObj.plotlist['Ub']

#Get profile
h, p = pObj.profile(axis=1, startrec=pObj.maxrec, endrec=pObj.maxrec)
h, Ub = UbObj.profile(axis=1, startrec=UbObj.maxrec, endrec=UbObj.maxrec)
p = p[:,0]
Ub = Ub[:,1]

#Sim box sizes
nx = pObj.Raw.ncx
ny = pObj.Raw.ncy
nz = pObj.Raw.ncz
xL = pObj.Raw.xL
yL = pObj.Raw.yL
zL = pObj.Raw.zL

# #Simulation parameters
with open('MD_dummy_fcc.py') as fObj:
	for l in fObj:
		if l.startswith('porousStart ='):
			porousStart = float(l.split()[-1])
		if l.startswith('porousEnd ='):
			porousEnd = float(l.split()[-1])
		if l.startswith('phi ='):
			phi = float(l.split()[-1])
		if l.startswith('dp ='):
			dp = float(l.split()[-1])
		if l.startswith('rho ='):
			rho = float(l.split()[-1])
		if l.startswith('mu ='):
			mu = float(l.split()[-1])
		if l.startswith('K ='):
			K = float(l.split()[-1])
		if l.startswith('Ubb ='):
			Ubb = float(l.split()[-1])
	
rp = 0.5*dp 
eps = 1 - phi
L = (porousEnd - porousStart - 1)*(yL/ny)
cvol = (xL/nx)*(yL/ny)*(zL/nz)
pvol = cvol*phi
cnp = pvol/((np.pi/6)*dp**3)
Cd = 3*np.pi*mu*dp*K
cCd = cnp*Cd*eps
pBot = 4.5*mu*phi*K*L*Ubb/(rp**2)

#Expected solution for velocity and pressure field
UbSol = np.concatenate(
    (Ubb*np.ones([int(porousStart), 1]), 
    (Ubb/eps)*np.ones([int(porousEnd)-int(porousStart), 1]), 
    Ubb*np.ones([ny-int(porousEnd), 1])))

pSol = np.concatenate(
    (pBot*np.ones([int(porousStart)]), 
    np.linspace(pBot, 0, int(porousEnd)-int(porousStart)), 
    np.zeros([ny-int(porousEnd)])))


def compare_results(x, xSol, tolrel=0.01, tolabs=1e-3):
	errx = np.zeros_like(xSol)
	for i in range(0, len(errx)):
		if xSol[i] == 0.:
			errx[i] = abs(xSol[i] - x[i]) <= tolabs
		else:
			errx[i] = abs((xSol[i] - x[i])/xSol[i]) <= tolrel

	return errx

errUb = compare_results(Ub, UbSol)
errp = compare_results(p, pSol)

#def test_velocityfield():
#	idx = np.where(errUb == False)[0]
#	print('Velocity profile at height h = ' + str(h[idx]) + ' exceeds the specified error tolerance.')
#	assert(errUb.all())

#def test_pressurefield():
#	idx = np.where(errp == False)[0]
#	print('Pressure profile at height h = ' + str(h[idx]) + ' exceeds the specified error tolerance.')
#	assert(errp.all())

data = np.loadtxt('regression_data/fcc_dummy_regression.txt')
UbReg = data[:,1]
pReg = data[:,2]
tolReg = 1e-3
errRegUb = abs(UbReg - Ub) <= tolReg
errRegp = abs(pReg - p) <= tolReg

def test_regression_velocityfield():
	idx = np.where(errRegUb == False)[0]
	print('Velocity profile at height h = ' + str(h[idx]) + ' fails the regression test.')
	assert(errRegUb.all())

def test_regression_pressurefield():
	idx = np.where(errRegp == False)[0]
	print('Pressure profile at height h = ' + str(h[idx]) + ' fails the regression test.')
	assert(errRegp.all())

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    #Plot velocity and pressure profile
    plt.plot(h, Ub, 'r-')
    plt.plot(h, UbSol, 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Velocity (cm/s)')
    plt.legend(('Numerical', 'Expected Solution'))
    plt.savefig('fig_velocityfield.png')
    plt.close()

    plt.plot(h, p, 'r-')
    plt.plot(h, pSol, 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Pressure (Ba = 0.1Pa)')
    plt.legend(('Numerical', 'Expected Solution'))
    plt.savefig('fig_pressurefield.png')
    plt.close()

except ImportError:
    pass
