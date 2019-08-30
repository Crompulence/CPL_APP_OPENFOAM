import os
import sys
import subprocess as sp
import numpy as np
import pytest

# Import postproclib
sys.path.insert(0, "./pyDataView/")
try:
    import postproclib as ppl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/pyDataView.git ./pyDataView"
    downloadout = sp.check_output(cmd, shell=True)
    sys.path.insert(0, "./pyDataView")
    import postproclib as ppl

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from LAMMPS_Input import LAMMPS_Input
from DragForce import DragForce, Stokes, DiFelice

def run_coupled(run_bash_script='run.sh'):
    try:
        cmd = './' + run_bash_script
        p = sp.Popen(cmd, 
            stdout=sp.PIPE, 
            stderr=sp.STDOUT, 
            shell=True)
        while p.poll() is None:
            l = p.stdout.readline()
            print(l.rstrip())
        print(p.stdout.read())
    except:
        print('Error running bash run script' + run_bash_script + ' in base directory.')
        raise

def read_print_data(file):
    data = np.loadtxt(file, skiprows=1)
    t = data[:,0]
    xy = data[:,1]

    return t, xy

def analytical_pressure(dragModel, muf, rhof, Uf, dp, epsf, g, xyzL, xyz_orig):

    # Obtain the pressure gradient across sample, and the inlet pressure
    # (assuming that the outlet pressure is zero)
    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=0., epsf=epsf)
    else:
        raise('Unknown drag force model specified')

    # Add pressure gradient contribution from gravity and calculate pressure
    # at the inlet
    gradP = fObj.calculate_pressure_gradient(epsf=epsf, Uf=Uf, Vp=0., dp=dp)
    pSol = (gradP + rhof*g)*(xyzL[1] - xyz_orig[1])
    
    return pSol
    
def get_velocity_pressure_porosity_profile(CFD_folder='./openfoam/', 
    velocity_var='Ub', pressure_var='p', porosity_var='alpha', axis_val=1):
    # Get plotting object
    PPObj = ppl.OpenFOAM_PostProc(CFD_folder)
    UObj = PPObj.plotlist[velocity_var]
    pObj = PPObj.plotlist[pressure_var]
    nObj = PPObj.plotlist[porosity_var]

    # Get profile (note that CPL calculates 1-eps)
    h, U = UObj.profile(axis=axis_val, startrec=UObj.maxrec, endrec=UObj.maxrec)
    h, p = pObj.profile(axis=axis_val, startrec=pObj.maxrec, endrec=pObj.maxrec)
    h, n = nObj.profile(axis=axis_val, startrec=nObj.maxrec, endrec=nObj.maxrec)
    h = h
    U = U[:,1]
    p = np.squeeze(p)
    n = np.squeeze(np.ones_like(n)-n)

    return h, U, p, n

def compare_pressure(pSol, p, h, xyz_orig, xyzL, tol=0.01):
    # Test pressure profile at the end of simulation.
    pSolProfile = np.interp(h, [xyz_orig[1],xyzL[1]], [pSol,0.])
    for i in range(len(pSolProfile)):
        err = (abs(pSolProfile[i] - p[i])/pSolProfile[i] <= tol)
        assert err, ('For h = {:.2f}, the obtained pressure of {:.6f} does not match analytical'.format(h[i], p[i])
                     + ' solution {:.6f} within {:.2f}% relative error'.format(pSolProfile[i], tol*100))

# ----- Main ----- #

# --------------------------------------------------------
# USER INPUT (ensure that these are the same as in MD_dummy_column.py) 
# --------------------------------------------------------
# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([0.1, 10.0, 0.1], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Parameters for model
dragModel = 'DiFelice'
muf = 1.e-2
rhof = 1.0
Uf = 0.1
epsf = 0.4764
Np_cell = 2 # Number of particles per fluid cell
dp = 0.1
g = 981.
# --------------------------------------------------------
# --------------------------------------------------------


@pytest.fixture(scope="module")
def setup():

    # Run coupled simulation
    run_coupled()

    # Extract input parameters from lammps input script
    pSol = analytical_pressure(dragModel, muf, rhof, Uf, dp, epsf, g, xyzL, xyz_orig)

    # Extract velocity, pressure and eps profile (at end of simulation only)
    h, U, p, n = get_velocity_pressure_porosity_profile()

    return pSol, h, U, p, n

tols = [1.0, 0.1, 0.01, 0.001, 0.0002]
@pytest.mark.parametrize("tols", tols)
def test_displacement(setup, tols):

    pSol, h, U, p, n = setup

    compare_pressure(pSol, p, h, xyz_orig, xyzL, tols)

if __name__ == "__main__":

    pSol, h, U, p, n = setup()

    import matplotlib.pyplot as plt

    # Plot Pressure Profile
    plt.plot(h, p, 'r-')
    plt.plot([xyz_orig[1],xyzL[1]], [pSol,0.], 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Pressure (0.1Pa)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    plt.savefig('fig_pressure.png')
    plt.close()
