import numpy as np
from mpi4py import MPI
import sys

# Import CPL library
from cplpy import CPL

# initialise MPI
comm = MPI.COMM_WORLD

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([2., 10., 2.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([10, 50, 10], order='F', dtype=np.int32)

# Initialise CPL
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

# Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=8, send_size=9)

# Main time loop
from StokesAnalytical import StokesAnalytical as SA
SAobj = SA()
for time in range(2501):

    print(time)
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    sendbuf[:,:,:,:] = 0
    sendbuf[4,:,:,:] = SAobj.rhof*SAobj.g

    CPL.send(sendbuf)

    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e] 
    recvbuf, ierr = CPL.recv(recvbuf)


CPL.finalize()
MPI.Finalize()

# Get numerical results and analytical solution
SAobj.get_analytical_solution()

# Test against analytical solution
SAobj.test_displacement()
SAobj.test_velocity()
SAobj.test_terminal()

# Plot results
SAobj.plot_displacement_velocity()