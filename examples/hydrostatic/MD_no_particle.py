import numpy as np
from mpi4py import MPI
import sys

#Import CPL library
from cplpy import CPL

#initialise MPI
comm = MPI.COMM_WORLD

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([2., 10., 2.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)

#Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

#Main time loop
for time in range(1):

    print(time)
    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)
    
    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e] 
    sendbuf[...] = 0
    sendbuf[7,:,:,:] = 0.
    CPL.send(sendbuf)

CPL.finalize()
MPI.Finalize()




