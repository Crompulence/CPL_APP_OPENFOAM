import numpy as np
from mpi4py import MPI
import sys

from cplpy import CPL

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)

# Parameters of the cpu topology (cartesian grid)
npxyz = [1, 1, 1]
NProcs = np.product(npxyz)
xyzL = np.array([16.795961913825074, 45.349097, 16.795961913825074], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

if (MD_COMM.Get_size() != NProcs):
    print("Non-coherent number of processes in MD ", MD_COMM.Get_size(),
            " no equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2])
    MPI.Abort(errorcode=1)

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_md(cart_comm, xyzL, xyz_orig)

#Setup send and recv buffers
recv_array, send_array = CPL.get_arrays(recv_size=3, send_size=4)

#Set velocity
U = 2.
N = 10

n = 0
for time in range(-2,199):

    # recv data 
    recv_array, ierr = CPL.recv(recv_array)

    # send data
    send_array[...] = 0.0
    send_array[0,:,:,:] = U*N #Sum of particles velocity (np.sin(2.*np.pi*time/10000.))*N
    send_array[3,:,:,:] = N    #Number of molecules

    CPL.send(send_array) 

    print("MDTime=", time)       
    
CPL.finalize()
MPI.Finalize()




