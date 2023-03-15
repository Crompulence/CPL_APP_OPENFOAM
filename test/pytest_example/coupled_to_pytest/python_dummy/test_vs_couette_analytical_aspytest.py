#!/usr/bin/env python
import numpy as np
import pytest

from CouetteAnalytical import CouetteAnalytical as CA

@pytest.fixture(scope="module")
def setup():

    #Import CPL library
    from cplpy import CPL

    #initialise MPI
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    #Check run as part of a coupled run
    comm.rank

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
    dt = 0.05
    U = 1.0
    nu = 1.004e-2
    Re = xyzL[1]/nu   #Note Reynolds in independent of velocity in analytical fn
    ncx = CPL.get("ncx")
    ncy = CPL.get("ncy")
    ncz = CPL.get("ncz")
    CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)

    #Yield statement delineates end of setup and start of teardown
    yield [CPL, MD_COMM, recvbuf, sendbuf, CAObj, dt, U, nu]
    CPL.finalize()
    MPI.Finalize()


#Main time loop
time = list(range(1000))
@pytest.mark.parametrize("time", time)
def test_loop(setup, time):

    #Get run paramenters from setup
    CPL, MD_COMM, recvbuf, sendbuf, CAObj, dt, U, nu = setup

    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e]
    sendbuf[...] = 0.
    CPL.send(sendbuf)

    #Get analytical solution
    y_anal, u_anal = CAObj.get_vprofile(time*dt)

    #Assert error bounds for L2 norm
    ur = np.mean(recvbuf[0,:,:,:],(0,2))
    error = np.sum(np.abs(100*(u_anal[1:-1:2] - ur)/U))
    print((time, "Error = ", error))
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





