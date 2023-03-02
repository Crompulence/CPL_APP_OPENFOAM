#!/usr/bin/env python
from mpi4py import MPI
from cplpy import CPL
import numpy as np

comm = MPI.COMM_WORLD
CPL = CPL()

MD_COMM = CPL.init(CPL.MD_REALM)

CPL.setup_md(MD_COMM.Create_cart([2, 1, 1]), 
             xyzL=[476.22031559046, 476.22031559046, 9.5244063118092], 
             xyz_orig=[0.0, 0.0, 0.0])

recv_array, send_array = CPL.get_arrays(recv_size=3, send_size=5)
uwall = 0.; vwall = 0.5
olap_limits = np.zeros(6); portion = np.zeros(6)
olap_limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(olap_limits)
dV = CPL.get("dx")*CPL.get("dy")*CPL.get("dz")
for time in range(501):
    recv_array, ierr = CPL.recv(recv_array)
    for i in range(send_array.shape[0]):
        for k in range(send_array.shape[2]):
            ig = i + portion[0]
            print(time, i, k, ig)
            if (ig > 60 and ig < 90):
                rho = 0.02;
                send_array[3,i,0,k] = rho*dV
                send_array[0,i,0,k] = uwall*send_array[3,i,0,k]
                send_array[1,i,0,k] = vwall*send_array[3,i,0,k]
            else:
                rho = 0.7;
                send_array[3,i,0,k] = rho*dV
                send_array[0,i,0,k] = uwall*send_array[3,i,0,k]
                send_array[1,i,0,k] = 0.0*send_array[3,i,0,k]
            send_array[4,i,0,k] = 0.95
    CPL.send(send_array)

    print("MD time", time)

CPL.finalize()
MPI.Finalize()


