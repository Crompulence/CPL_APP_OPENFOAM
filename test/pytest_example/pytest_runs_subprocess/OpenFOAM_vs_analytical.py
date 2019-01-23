import numpy as np
import sys
import subprocess as sp

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


def test_error(error, time):
    if time < 10:
        assert error < 10., "Error in inital 10 steps greater than 10%"
    elif time < 30:
        assert error < 5., "Error between 10 and 30 steps greater than 5%"
    elif time < 50:
        assert error < 3., "Error between 30 and 50 steps greater than 3%"
    elif time < 300:
        assert error < 2., "Error between 50 and 300 steps greater than 2%"
    elif time < 500:
        assert error < 0.2, "Error between 300 and 500 steps greater than 0.2%"
    else:
        assert error < 0.1, "Error after 500 steps greater than 0.1%"


def check_OpenFOAM_vs_Analytical(fdir, plotstuff = False):

    OpenFOAMfdir = fdir+"/cfd_data/openfoam/"
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
    with open(scriptdir+"test_vs_couette_analytical.py", 'r') as f:
        for l in f:
            if "U =" in l:
                d = l.split("=")
                U = float(d[1].replace("\n","").replace(" ",""))
    nu = OpenFOAMuObj.Raw.nu[0]
    Re = (xyzL[1]+dy/2.)/nu
    CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1)

    if plotstuff:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
    n = 0
    for time in range(enditer):

        #Plot data
        if time%OpenFOAMwriteinterval == 0:
            rec = int(time/float(OpenFOAMwriteinterval))
            try:
                OpenFOAMpp = ppl.OpenFOAM_PostProc(OpenFOAMfdir)
                OpenFOAMuObj = OpenFOAMpp.plotlist['Ub']#ppl.OpenFOAM_vField(OpenFOAMfdir, parallel_run=True)
                y, u = OpenFOAMuObj.profile(1,startrec=rec,endrec=rec)
                y_anal, u_anal = CAObj.get_vprofile(time*dt, flip=False)
                error = (u_anal[1:-1:2] - u[:,0])/u_anal[1:-1:2]

                if plotstuff:
                    l, = ax.plot(u[:,0], y, 'ro-', label="OpenFOAM domain from file")
                    ax.plot(u_anal[1:-1:2],y_anal[1:-1:2], 'k.-', label="Analytical Solution")

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
                    print("Time = ", time, "Error= ", l2)
                    test_error(l2, time)

            except AssertionError as e:
                print("AssertionError ", e)
                raise

            except:# ppl.field.OutsideRecRange:
                print("Error result missing", time, OpenFOAMuObj.maxrec, rec)
                raise


if __name__ == "__main__":
    check_OpenFOAM_vs_Analytical("./run1.0/", plotstuff=True)
