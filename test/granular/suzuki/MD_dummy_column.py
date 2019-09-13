import sys
import numpy as np
from mpi4py import MPI
from cplpy import CPL

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from DragForce import DragForce, Stokes, DiFelice

# --------------------------------------------------------
# USER INPUT 
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
Np_cell = 1. # Number of particles per fluid cell
dp = 0.1
g = 981.
# --------------------------------------------------------
# --------------------------------------------------------
print('initialise')

# initialise MPI
comm = MPI.COMM_WORLD

# Initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)

# Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

# Setup drag force object
print('drag')
if dragModel == 'Drag' or dragModel == 'Stokes':
    fObj = Stokes(muf=muf, dp=dp)
elif dragModel == 'DiFelice':
    fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=0., epsf=epsf)
else:
    raise('Unknown drag force model specified')

# Main time loop
print('time loop')
for time in range(200):

    print(time)
    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    # Zero send buffer. Set the sum of drag coeff and particle volumes only
    # (noting that CPLSocketFOAM divides by sample volume). Assumes entire
    # column is filled with particles.
    # [UxSum, UySum, UzSum, FxSum, FySum, FzSum, CdSum, VolSum] 
    sendbuf[:,:,:,:] = 0.
    sendbuf[6,:,:,:] = Np_cell*fObj.B
    sendbuf[7,:,:,:] = Np_cell*(np.pi/6)*(dp**3)

    CPL.send(sendbuf)

CPL.finalize()
MPI.Finalize()