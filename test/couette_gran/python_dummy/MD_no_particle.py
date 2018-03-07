import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import sys


from cplpy import CPL
from CouetteAnalytical import CouetteAnalytical as CA

import sys
sys.path.append("/home/es205/codes/python/pyDataView")
import postproclib as ppl

#Import CPL library
from cplpy import CPL

OpenFOAMfdir = "./openfoam"

#Define simulatiom parameters
dt = 1e-5
g = 9.81
D = 1.e-4
rho_Gran = 2.65e+03
mass = rho_Gran*(4./3.)*np.pi*(D/2.)**3
mu = 1.004e-2
timestep_ratio = 100

dt = 0.1
tplot = 1.
OpenFOAMwriteinterval = tplot/dt


#initialise MPI
comm = MPI.COMM_WORLD
if (comm.Get_size() != 1):
    print("Simplified code only works for one proccesor!")
    MPI.Abort(errorcode=1)

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([2., 2., 2.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)


#initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
#CPL.set_timing(0, 0, dt)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)

#Setup send and recv buffers
cnst_limits = CPL.get_cnst_limits()
cnst_portion = CPL.my_proc_portion(cnst_limits)
[cnst_ncxl, cnst_ncyl, cnst_nczl] = CPL.get_no_cells(cnst_portion)
recvbuf = np.zeros((9, cnst_ncxl, cnst_ncyl, cnst_nczl), order='F', dtype=np.float64)

olap_limits = CPL.get_olap_limits()
BC_limits = np.array([olap_limits[0],olap_limits[1],
                      olap_limits[2],olap_limits[3],
                      olap_limits[4],olap_limits[5]], order='F', dtype=np.int32)
BC_portion = CPL.my_proc_portion(BC_limits)
[BC_ncxl, BC_ncyl, BC_nczl] = CPL.get_no_cells(BC_portion)

#Plot output
fig, ax = plt.subplots(1,1)
plt.ion()
plt.show()

#Analytical solution
U = 1.
nu = 1e-02
Re = U/nu
ncy = CPL.get("ncy")
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax = xyzL[1], npoints=ncy)

#Porous region solution
pcell = 7
dy = CPL.get("dy")
yl_cfd = CPL.get("yl_cfd")

CApObj = CA(Re=Re, U=U, Lmin=yl_cfd-(pcell+2.5)*dy, Lmax = xyzL[1], npoints=ncy)
#y_p = np.linspace(yl_cfd-(pcell+2)*dy,yl_cfd, ncy)

# Coupling parameters
eps = 0.0001
A = 1./eps
l2 = []

#Define a field of random eps
eps = 0.4 + 0.1*np.random.randn(cnst_ncxl, cnst_ncyl, cnst_nczl)

#Main time loop
pcount = 0
for time in range(1000):

    print(time)
    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf, cnst_portion)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e]
    sendbuf = np.zeros([8, cnst_ncxl, cnst_ncyl, cnst_nczl], order='F', dtype=np.float64)
    sendbuf[7,:,:,:] = 1.0

    #Set pcell to zero porosity
    #sendbuf[7,:,pcell:pcell+2,:] = eps
    #sendbuf[7,:,pcell+2:pcell+4,:] = 1.0

    # Set porous force based on drag difference between velocity
    #sendbuf[3,:,pcell+2:pcell+4,:] = A*( sendbuf[0,:,pcell+2:pcell+4,:]
    #                                    -recvbuf[0,:,pcell+2:pcell+4,:]) 

#    dedy = np.gradient(sendbuf[7,:,:,:],axis=1)
#    dudy = np.gradient(recvbuf[0,:,:,:],axis=1)
#    F = mu*dedy*dudy/(eps*dy**2)
#    print(dedy[1,pcell:pcell+4,1], dudy[1,pcell:pcell+4,1], F[1,pcell:pcell+4,1])
#    sendbuf[3,:,pcell+2:pcell+4,:] = F[:,pcell:pcell+1,:]
    CPL.send(sendbuf, BC_limits)

    #Plot data
    if time%OpenFOAMwriteinterval == 0:

        rec = int(time/OpenFOAMwriteinterval)

        try:

            #Plot values from OpenFOAM read by postproc library 
            OpenFOAMuObj = ppl.OpenFOAM_vField(OpenFOAMfdir)
            y, u = OpenFOAMuObj.profile(1,startrec=rec,endrec=rec)
            l, = ax.plot(u[:,0], y, 'ro-', label="OpenFOAM domain from file")

            #halou = OpenFOAMuObj.read_halo(startrec=rec,endrec=rec)
            #ax.plot(np.mean(halou[:,:,:,:,0],(0,2)), y[0]-0.5*dy, 'bs',ms=10, label="OpenFOAM halo from file")

            #Plot values recieved from SediFOAM
            ur = np.mean(recvbuf[0,:,:,:],(0,2))
            ax.plot(ur, y,'g^',ms=10,label="Recieved value from SediFOAM")

            #Plot analytical solution
            y_anal, u_anal = CAObj.get_vprofile(time*dt)
            ax.plot(u_anal[:],y_anal, 'k.-', label="Analytical Solution")

            y_p, u_p = CApObj.get_vprofile(time*dt)
            ax.plot(u_p, y_p, 'kx--', label="Analytical Solution porous")


            error = (u_p[::2] - ur[pcell+3:])/u_p[::2]
            l2.append(np.sqrt(np.sum(error[1:]**2)))

            ax.set_xlim([-0.1,1.1])
            ax.set_xlabel("$u$",fontsize=24)
            ax.set_ylabel("$y$",fontsize=24)

            plt.legend(loc=4)
            plt.pause(0.001)

         
            #plt.savefig('out{:05}.png'.format(pcount));
            pcount += 1

            ax.cla()
            print("MDTime=", time, OpenFOAMwriteinterval, OpenFOAMuObj.maxrec, rec)
        except: #ppl.field.OutsideRecRange:
            #print("Error result missing", time, OpenFOAMuObj.maxrec, rec)
            raise


print(l2)
#import cPickle as pickle
#pickle.dump(l2, open("eps"+str(eps)+".p", "w+"))

CPL.finalize()
MPI.Finalize()




