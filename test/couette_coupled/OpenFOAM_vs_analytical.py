import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/es205/codes/python/pyDataView")
import postproclib as ppl

sys.path.append("/home/es205/codes/cpl/cpl-library/utils/")
from CouetteAnalytical import CouetteAnalytical as CA

def test_error(error, time):
    if time < 350:
        assert error < 8., "Error in inital 350 steps greater than 8%"
    elif time < 600:
        assert error < 1., "Error between 350 and 600 steps greater than 1%"
    elif time < 1000:
        assert error < 0.1, "Error after 1000 steps greater than 0.1%"

def check_OpenFOAM_vs_Analytical(fdir, plotstuff = False):

    OpenFOAMfdir = fdir+"/cfd_data/openfoam_ico/"
    print("Openfoam dir = ", OpenFOAMfdir)
    OpenFOAMuObj = ppl.OpenFOAM_vField(OpenFOAMfdir, parallel_run=True)
    dt = float(OpenFOAMuObj.Raw.delta_t)
    tplot = float(OpenFOAMuObj.Raw.header.headerDict['controlDict']['writeInterval'])
    endtime = float(OpenFOAMuObj.Raw.header.headerDict['controlDict']['endTime'])
    enditer = int(endtime/dt)
    OpenFOAMwriteinterval = int(tplot/dt)

    # Parameters of the cpu topology (cartesian grid)
    xyzL = [OpenFOAMuObj.Raw.xL, OpenFOAMuObj.Raw.yL, OpenFOAMuObj.Raw.zL]
    dx = OpenFOAMuObj.Raw.dx
    dy = OpenFOAMuObj.Raw.dy
    dz = OpenFOAMuObj.Raw.dz

    ncx = OpenFOAMuObj.Raw.ncx
    ncy = OpenFOAMuObj.Raw.ncy
    ncz = OpenFOAMuObj.Raw.ncz

    #Analytical solution
    scriptdir = fdir+"/md_data/python_dummy/"
    with open(scriptdir+"MD_sendrecv.py", 'r') as f:
        for l in f:
            if "U =" in l:
                d = l.split("=")
                U = float(d[1].replace("\n","").replace(" ",""))
    nu = OpenFOAMuObj.Raw.nu[0]
    Re = (xyzL[1]+dy/2.)/nu
    CAObj = CA(Re=Re, U=U, Lmin=0.-dy/2., Lmax=xyzL[1], npoints=2*ncy+2)

    if plotstuff:
        fig, ax = plt.subplots(1,1)
    n = 0
    for time in range(enditer):

        #Plot data
        if time%OpenFOAMwriteinterval == 0:
            rec = int(time/float(OpenFOAMwriteinterval))
            try:
                OpenFOAMuObj = ppl.OpenFOAM_vField(OpenFOAMfdir, parallel_run=True)
                y, u = OpenFOAMuObj.profile(1,startrec=rec,endrec=rec)
                halou = OpenFOAMuObj.read_halo(startrec=rec,endrec=rec, haloname="CPLReceiveMD")
                halout = OpenFOAMuObj.read_halo(startrec=rec,endrec=rec, haloname="movingWall")
                y_anal, u_anal = CAObj.get_vprofile(time*dt, flip=True)
                error = (u_anal[2:-1:2] - u[:,0])/u_anal[2:-1:2]

                if plotstuff:
                    l, = ax.plot(u[:,0], y, 'ro-', label="OpenFOAM domain from file")
                    ax.plot(np.mean(halou[:,:,:,:,0],(0,2)), y[0]-0.5*dy, 'bs',ms=10, label="OpenFOAM halo from file")
                    ax.plot(np.mean(halout[:,:,:,:,0],(0,2)), y[-1]+0.5*dy, 'bs',ms=10, label="OpenFOAM halo fixed")
                    ax.plot(u_anal,y_anal, 'k.-', label="Analytical Solution")

                    #ax.plot(10.*(u[:,0]-u_anal[-2:0:-2]),y,'y--')
                    ax.set_xlim([-0.1,U*1.1])
                    ax.set_xlabel("$u$",fontsize=24)
                    ax.set_ylabel("$y$",fontsize=24)

                    plt.legend(loc=1)
                    plt.pause(0.001)
                    plt.savefig('out{:05}.png'.format(n)); n+=1

                    ax.cla()

            
                l2 = np.sqrt(np.sum(error[1:]**2))
                if not np.isnan(l2):
                    print(time, l2)
                    test_error(l2, time)

            except AssertionError as e:
                print("AssertionError ", e)
                raise

            except:# ppl.field.OutsideRecRange:
                print("Error result missing", time, OpenFOAMuObj.maxrec, rec)
                raise

