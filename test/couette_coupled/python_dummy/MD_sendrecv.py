import numpy as np
from mpi4py import MPI

from cplpy import CPL
from draw_grid import draw_grid

import sys
sys.path.append("/home/es205/codes/python/pyDataView")
import postproclib as ppl

sys.path.append("/home/es205/codes/cpl/cpl-library/utils/")
from CouetteAnalytical import CouetteAnalytical as CA

OpenFOAMfdir = "./openfoam_ico/"
dt = 0.005
tplot = 0.16
OpenFOAMwriteinterval = tplot/dt

N = 10

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
rank = MD_COMM.Get_rank()
nprocs_realm = MD_COMM.Get_size()

# Parameters of the cpu topology (cartesian grid)
npxyz = [1, 1, 1]
NProcs = np.product(npxyz)
xyzL = np.array([195.2503206, 18.62550553, 133.3416884], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

if (nprocs_realm != NProcs):
    print("Non-coherent number of processes in MD ", nprocs_realm,
            " no equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2])
    MPI.Abort(errorcode=1)

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_md(cart_comm, xyzL, xyz_orig)
ir, jr, kr = cart_comm.Get_coords(cart_comm.Get_rank())
ncx_l = CPL.get("ncx")/npxyz[0]
ncy_l = CPL.get("ncy")/npxyz[1]
ncz_l = CPL.get("ncz")/npxyz[2]

# === Plot both grids ===
dx = CPL.get("xl_md")/float(CPL.get("ncx"))
dy = CPL.get("yl_md")/float(CPL.get("ncy"))
dz = CPL.get("zl_md")/float(CPL.get("ncz"))

#Analytical solution
U = 2.
nu = 10.
Re = (xyzL[1]+dy/2.)/nu
ncy = CPL.get("ncy")
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1)

#Setup send and recv buffers
cnst_limits = CPL.get_cnst_limits();
cnst_portion = CPL.my_proc_portion(cnst_limits)
[cnst_ncxl, cnst_ncyl, cnst_nczl] = CPL.get_no_cells(cnst_portion)
recv_array = np.zeros((3, cnst_ncxl, cnst_ncyl, cnst_nczl), order='F', dtype=np.float64)

olap_limits = CPL.get_olap_limits()
BC_limits = np.array(CPL.get_bnry_limits(), dtype=np.int32)
BC_portion = CPL.my_proc_portion(BC_limits)
[BC_ncxl, BC_ncyl, BC_nczl] = CPL.get_no_cells(BC_portion)

n = 0
for time in range(-2,1999):

    # recv data 
    recv_array, ierr = CPL.recv(recv_array, cnst_portion)

    # send data
    send_array = np.zeros((5, BC_ncxl, BC_ncyl, BC_nczl), order='F', dtype=np.float64)
    send_array[0,:,:,:] = U*N #Sum of particles velocity (np.sin(2.*np.pi*time/10000.))*N
    send_array[3,:,:,:] = N    #Number of molecules

    CPL.send(send_array, BC_portion) 

    print("MDTime=", time)       
    
CPL.finalize()
MPI.Finalize()




