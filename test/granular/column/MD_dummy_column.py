import sys
import numpy as np
from mpi4py import MPI
from cplpy import CPL

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from OpenFOAM_Input import OpenFOAM_Input
from DragForce import DragForce, Stokes, DiFelice, Ergun

# --------------------------------------------------------
# USER INPUT - START
# --------------------------------------------------------
npxyz = [1, 1, 1]
xyzL = [0.1, 10.0, 0.1]
xyz_orig = [0.0, 0.0, 0.0]

dragModel = 'Ergun'
epsf = 0.476401224402
Np_cell = 1.
dp = 0.1
# --------------------------------------------------------
# USER INPUT - END
# --------------------------------------------------------
cObj = OpenFOAM_Input(case_dir='./openfoam/')
_, Uf = cObj.read_0_file('./openfoam/0/Ub', 'inlet')
rhof = cObj.read_constant_file('./openfoam/constant/transportProperties', 'rhob')
nuf = cObj.read_constant_file('./openfoam/constant/transportProperties', 'nub')
Uf = Uf[1]
muf = nuf*rhof

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array(npxyz, order='F', dtype=np.int32)
xyzL = np.array(xyzL, order='F', dtype=np.float64)
xyz_orig = np.array(xyz_orig, order='F', dtype=np.float64)

# initialise MPI
comm = MPI.COMM_WORLD

# Initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)

# Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

# Setup drag force object
if dragModel == 'Drag' or dragModel == 'Stokes':
    fObj = Stokes(muf=muf, epsf=epsf, dp=dp)
elif dragModel == 'DiFelice':
    fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
elif dragModel == 'Ergun':
    fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
else:
    raise ValueError('Unknown drag force model specified')

# Main time loop
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
    sendbuf[6,:,:,:] = Np_cell*fObj.B*epsf
    sendbuf[7,:,:,:] = Np_cell*(np.pi/6)*(dp**3)

    CPL.send(sendbuf)

CPL.finalize()
MPI.Finalize()