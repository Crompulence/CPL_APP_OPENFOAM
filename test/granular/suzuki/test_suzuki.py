import os
import sys
import errno
import pytest
import subprocess as sp
import numpy as np

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
from MOCK_Input import MOCK_Input, MOCK_Writer
from OpenFOAM_Input import OpenFOAM_Input, OpenFOAM_Writer_0_File
from DragForce import DragForce, Stokes, DiFelice, Ergun

# Run coupled simulation as subprocess
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
        raise RuntimeError('Error running bash run script' + run_bash_script + ' in base directory.')

# Extract input parameters from MD/DEM mock script and OpenFOAM case. Check
# that these are consistent.
def get_input_parameters(md_input_file='MD_dummy_suzuki.py', cfd_case_dir='./openfoam/'):
    # MD/DEM mock script
    mObj = MOCK_Input(md_input_file)

    # OpenFOAM case 
    cObj = OpenFOAM_Input(cfd_case_dir)
    _, Uf = cObj.read_0_file('./openfoam/0/Ub', 'inlet')
    nuf = cObj.read_constant_file('./openfoam/constant/transportProperties', 'nub')
    mObj.rhof = cObj.read_constant_file('./openfoam/constant/transportProperties', 'rhob')
    mObj.g = cObj.read_constant_file('./openfoam/constant/environmentalProperties', 'g')[1]
    mObj.Uf = Uf[1]
    mObj.muf = nuf*mObj.rhof
    
    return mObj

# Set the input parameters for the simulations. At present, only the inlet
# fluid velocity and the drag model can be adjusted.
def set_input_parameters(Uf, dragModel, 
        md_input_file='./MD_dummy_suzuki.py', cfd_case_dir='./openfoam/'):
    
    MOCK_Writer(md_input_file, 'dragModel', dragModel)
    OpenFOAM_Writer_0_File(cfd_case_dir + '0/Ub', 'inlet', Uf, isScalar=False, axis_val=1)

# Calculate the analytical solution for the pressure at the inlet, assuming
# that the pressure at the outlet is zero.
def analytical_pressure(mObj):
    
    dragModel = mObj.dragModel
    muf = mObj.muf
    rhof = mObj.rhof
    Uf = mObj.Uf
    dp = mObj.dp
    epsf = mObj.epsf
    xyzL = mObj.xyzL
    xyz_orig = mObj.xyz_orig

    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, epsf=epsf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    else:
        raise ValueError('Unknown drag force model specified')

    gradP = fObj.calculate_pressure_gradient(Uf=Uf/epsf, Vp=0.)
    pSol = gradP*(xyzL[1] - xyz_orig[1])
    
    return pSol

# Obtain the pressure profile from the numerical simulation
def get_pressure_profile(CFD_folder='./openfoam/', pressure_var='p', axis_val=1):
    PPObj = ppl.OpenFOAM_PostProc(CFD_folder)
    pObj = PPObj.plotlist[pressure_var]
    h, p = pObj.profile(axis=axis_val, startrec=pObj.maxrec, endrec=pObj.maxrec)
    p = np.squeeze(p)

    return h, p

# Plot pressure profile obtained from numerical simulation and analytical
# solution. Save the file in the results directory (which is created if
# required) and also save the data used for plotting as .npz file.
def plot_pressure(h, p, pSol, xyz_orig, xyzL, file_name='fig_pressure'):
    plt.plot(h, p, 'r-')
    plt.plot([xyz_orig[1],xyzL[1]], [pSol,0.], 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Pressure (0.1Pa)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    plt.savefig(file_name + '.png')
    np.savez(file_name + '.npz', h=h, p=p, pSol=pSol, xyz_orig=xyz_orig, xyzL=xyzL)
    plt.close()

# Compare the pressure along the length of the sample at the end of the
# simulation for a specified relative error.
def compare_pressure(h, p, pSol, xyz_orig, xyzL, tol=0.01):
    pSolProfile = np.interp(h, [xyz_orig[1],xyzL[1]], [pSol,0.])
    for i in range(len(pSolProfile)):
        err = (abs(pSolProfile[i] - p[i])/pSolProfile[i] <= tol)
        assert err, ('For h = {:.2f}, the obtained pressure of {:.6f} does not match analytical'.format(h[i], p[i])
                     + ' solution {:.6f} within {:.2f}% relative error'.format(pSolProfile[i], tol*100))

# ----- Main ----- #
dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.1, 0.2, 0.3, 0.4, 0.5]
@pytest.mark.parametrize('dragModel', dragModels)
@pytest.mark.parametrize('Uf', Uf_values)
def test_pressure(Uf, dragModel, plot_results=False):
    # Set input parametesr
    set_input_parameters(Uf, dragModel)

    # Run coupled simulation
    run_coupled()

    # Extract input parameters
    mObj = get_input_parameters()

    # Calculate analytical pressure profile
    pSol = analytical_pressure(mObj)

    # Extract numerical pressure profile (end of simulation)
    h, p = get_pressure_profile()

    # Correct for the presence of a buoyancy gradient
    p = p + mObj.rhof*mObj.g*np.flipud(h)

    # Plot the results
    if plot_results:
        import matplotlib.pyplot as plt
        plot_pressure(h, p, pSol, mObj.xyz_orig, mObj.xyzL, 
            file_name='./results/fig_pressure_Uf_{}_{}'.format(Uf, dragModel))

    compare_pressure(h, p, pSol, mObj.xyz_orig, mObj.xyzL, tol=0.02)

if __name__ == "__main__":
    test_pressure(Uf=0.5, dragModel='DiFelice')

    

