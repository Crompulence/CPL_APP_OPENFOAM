import numpy as np
import sys
import subprocess as sp
import time

# Import symwraplib
sys.path.insert(0, "./pyDataView/")
try:
    import postproclib as ppl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/pyDataView.git ./pyDataView"
    downloadout = sp.check_output(cmd, shell=True)
    sys.path.insert(0, "./pyDataView")
    import postproclib as ppl

from CouetteAnalytical import CouetteAnalytical as CA

def test_error(error, t):
    if t < 35:
        assert error < 0.03, "Error in inital 30 steps greater than 30%"
    elif t < 60:
        assert error < 0.01, "Error between 30 and 60 steps greater than 1%"
    elif t > 60:
        assert error < 0.008, "Error after 1000 steps greater than 1.0%"

def check_OpenFOAM_vs_Analytical(fdir, plotstuff = False, parallel_run=True):

    OpenFOAMfdir = fdir+"/cfd_data/openfoam_ico/"
    print("Openfoam dir = ", OpenFOAMfdir)
    OpenFOAMuObj = ppl.OpenFOAM_vField(OpenFOAMfdir, parallel_run=parallel_run)
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
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        if "dynamic" in plotstuff or plotstuff == True:
            plotstuff = "dynamic"
            recds = list(range(enditer))
        elif "summary" in plotstuff:
            recds = [5, int(0.1*enditer), int(0.25*enditer), enditer-1]
    else:
        plotstuff = "False"

    n = 0
    for t in range(enditer):

        #Plot data
        if t%OpenFOAMwriteinterval == 0:
            rec = int(t/float(OpenFOAMwriteinterval))
            try:
                read_attempt=5
                while True:
                    try:
                        OpenFOAMuObj = ppl.OpenFOAM_vField(OpenFOAMfdir, parallel_run=parallel_run)
                        y, u = OpenFOAMuObj.profile(1,startrec=rec,endrec=rec)
                        halou = OpenFOAMuObj.read_halo(startrec=rec,endrec=rec, haloname="CPLReceiveMD")
                        print(halou.shape, halou[...,0].min(),halou[...,0].max())
                        halout = OpenFOAMuObj.read_halo(startrec=rec,endrec=rec, haloname="movingWall")
                        break
                    except ppl.field.OutsideRecRange:
                        if (read_attempt == 0):
                            raise
                        else:
                            print("At record=", rec, "data not found", 
                                  "read attempts left=", read_attempt, 
                                  ". Waiting 2 seconds and trying again.")
                            time.sleep(2.)
                            read_attempt =- 1
                            continue
                            
                y_anal, u_anal = CAObj.get_vprofile(t*dt, flip=True)
                error = (u_anal[2:-1:2] - u[:,0])/U #/u_anal[2:-1:2]
                #error[u_anal[2:-1:2] < 0.005] = 0.
                if plotstuff != "False":
                    if t in recds:
                        l, = ax.plot(u[:,0], y, 'ro-', 
                                     label="OpenFOAM domain from file")
                        ax.plot(np.mean(halou[:,:,:,:,0],(0,2)), y[0]-0.5*dy, 
                                        'bs',ms=10, label="OpenFOAM halo from file")
                        ax.plot(np.mean(halout[:,:,:,:,0],(0,2)), y[-1]+0.5*dy, 
                                        'bs',ms=10, label="OpenFOAM halo fixed")
                        ax.plot(u_anal,y_anal, 'k.-', label="Analytical Solution")
                        ax.plot(-error,y_anal[2:-1:2], 'g--', label="Analytical Solution")

                    #ax.plot(10.*(u[:,0]-u_anal[-2:0:-2]),y,'y--')
                    if "dynamic" in plotstuff:
                        ax.set_xlim([-0.1,U*1.1])
                        ax.set_xlabel("$u$",fontsize=24)
                        ax.set_ylabel("$y$",fontsize=24)

                        plt.legend(loc=1)
                        plt.pause(0.001)
                        plt.savefig('out{:05}.png'.format(n)); n+=1

                        ax.cla()
            
                l2 = np.sqrt(np.sum(error[1:]**2))
                if not np.isnan(l2):
                    print(t, "L2 norm error =" , l2)
                    test_error(l2, t)

            except AssertionError as e:
                print("AssertionError ", e)
                raise

            except ppl.field.OutsideRecRange:
                print("Error result missing", t, OpenFOAMuObj.maxrec, rec)
                raise

    if "summary" in plotstuff:
        ax.set_xlim([-0.1,U*1.1])
        ax.set_xlabel("$u$",fontsize=24)
        ax.set_ylabel("$y$",fontsize=24)

        #plt.legend(loc=1)
        plt.savefig('summary.png')


if __name__ == "__main__":
    check_OpenFOAM_vs_Analytical("./run1.0/", plotstuff="dynamic")
