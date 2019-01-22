
# Subprocess Run Tests

This directory is an example for a base script for running a range of test cases using Pytest

## What does it do?

The range of cases are run include a number of different solvers and a number of different wall sliding speeds.
This is parameterised using pytest in test_couette.py as follows,

```python
#Parameterise range of cases
params = []
Uwall = [0.2, 0.5, 1.0, 2.0]
for exe in ["CPLSediFOAM", "CPLCFDDEMFoam"]:
    for u in Uwall:
        params.append([u, exe])
@pytest.mark.parametrize("wallvel, executable", params)
def test_newtest(wallvel, executable):

```
where the actual changes (wallvel) to the CFD input and the mock script are
made using simwraplib, which creates directories for each run and edits
 key input values in the mock scripts and OpenFOAM (see simwraplib 
[website][https://github.com/edwardsmith999/SimWrapPy] for more details).
The actual coupled run uses both openFOAM and mock script, 
where the mock Python script checks the OpenFOAM output against 
the analytical solution for Couette flow, as follows,

```python
import numpy as np
from mpi4py import MPI
import sys

from CouetteAnalytical import CouetteAnalytical as CA

#Import CPL library
from cplpy import CPL

#initialise MPI
comm = MPI.COMM_WORLD

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([1., 1., 1.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

#Analytical solution
plot = False
error_check=True
dt = 0.05
U = 1.
nu = 1.004e-2
Re = xyzL[1]/nu   #Note Reynolds in independent of velocity for some reason
ncx = CPL.get("ncx")
ncy = CPL.get("ncy")
ncz = CPL.get("ncz")
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)
dx = xyzL[0]/float(ncx)
dy = xyzL[1]/float(ncy)
dz = xyzL[2]/float(ncz)
dV = dx*dy*dz
y = np.linspace(dy/2., xyzL[1] - dy/2., num=ncy)

#Main time loop
pcount = 0
for time in range(1000):

    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e]
    sendbuf[...] = 0.
    CPL.send(sendbuf)

    ur = np.mean(recvbuf[0,:,:,:],(0,2))
    y_anal, u_anal = CAObj.get_vprofile(time*dt)

    #Plot values recieved from SediFOAM
    if error_check:
        error = np.sum(np.abs(100*(u_anal[1:-1:2] - ur)/U))
        print(time, "Error = ", error)
        try:
            test_error(error)
        except AssertionError as e:
            print("AssertionError ", e)
            comm.Abort()

CPL.finalize()
MPI.Finalize()
```

Any assertion error, raised by test_error, will stop the code and this is caught
 by the test_couette.py code and raised as an exception which is counted as
an error by pytest.

## What is tested?

This checks that a range of OpenFOAM solvers reproduce time evolving Couette flow for a range of wall speeds.
Couette flow is useful as it is one of the cases we have an analytical solutions for, but the actual checks to
the solver itself is limited to just the diffusion term.

## If it breaks, what do you check

There is a known problem with the interplay between python subprocesses and catching error codes (especially on Travis). 
It may be worth running this again to ensure tests are working as expected.
The error bounds in test_vs_couette_analytical.py under test_error are hardwired based on experience with the current 
solvers, 

```python
def test_error(error):
    if time < 10:
        assert error < 20., "Error in inital 10 steps greater than 20%"
    elif time < 30:
        assert error < 10., "Error between 10 and 30 steps greater than 10%"
    elif time < 50:
        assert error < 5., "Error between 30 and 50 steps greater than 5%"
    elif time < 300:
        assert error < 3., "Error between 50 and 300 steps greater than 3%"
    elif time < 500:
        assert error < 2., "Error between 300 and 500 steps greater than 2%"
    else:
        assert error < 1., "Error after 500 steps greater than 1%"

```

If you plan to use for another case, these will almost certainly need to be adapted to acceptable bounds for your problem.
