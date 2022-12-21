#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

from mpi4py import MPI
from cplpy import CPL

#Import for loading data and plotting
try:
    ppdir = '/home/es205/codes/python/pyDataView/'
    sys.path.append(ppdir)
    import postproclib as ppl
except ImportError:
    import subprocess as sp
    cmd = "git clone https://github.com/edwardsmith999/pyDataView.git ./pyDataView"
    downloadout = sp.check_output(cmd, shell=True)
    print(downloadout)
    sys.path.insert(0, "./pyDataView")
    import postproclib as ppl

if __name__ == '__main__':

    plot = False
    plot_CFD = False
    plot_vmd = False
    plot_grid = False
    plot_quiver = True
    savefig = False

    if plot_grid:
        from draw_grid import draw_grid

    #========================================================================
    # Load MD data
    #Get Post Proc Object
    fdir = '/home/es205/codes/flowmol/runs/Boiling_Square_Notch/Young/results/'
    startrec = 750; endrec = 1600
    yloc = 32
    MDPPObj = ppl.MD_PostProc(fdir)

    #Get plotting object
    rhoObj = MDPPObj.plotlist['rho']
    uObj = MDPPObj.plotlist['u']
    TObj = MDPPObj.plotlist['T']
    rho = rhoObj.read(startrec=startrec, endrec=startrec+1)
    Lx_MD = float(rhoObj.header.globaldomain1)
    Ly_MD = float(rhoObj.header.globaldomain2)
    Lz_MD = float(rhoObj.header.globaldomain3)

    #dy =float(rhoObj.header.binsize2)
    x = np.linspace(0.,Lx_MD,rho.shape[0])
    y = np.linspace(0.,Ly_MD,rho.shape[1])
    X,Y = np.meshgrid(x,y)

    dx_MD =float(rhoObj.header.binsize1)
    dy_MD =float(rhoObj.header.binsize2)

    #Now we have a psf file, better to use this (access to bond, etc)
    if plot and plot_vmd:
        try:
            import MDAnalysis
            mols = MDAnalysis.Universe(fdir+'polymer_topol.psf', fdir+'vmd_out.dcd')
            DCDObj = mols.trajectory
        except ImportError:
            print("MDAnalysis not available, no molecular data can be plotted")
            plot_vmd = False
        except IOError:
            print("polymer_topol.psf not found, trying vmd_out.dcd")
            try:
                DCDObj = MDAnalysis.coordinates.DCD.DCDReader(fdir+"vmd_out.dcd")
            except IOError:
                print("vmd_out.dcd not available,  no molecular data can be plotted")

    #CFD variables and load data to plot
    if plot and plot_CFD:
        CFDPPObj = ppl.All_PostProc("./")
        print(CFDPPObj)
        CFDrhoObj = CFDPPObj.plotlist['rho']
        CFDUObj = CFDPPObj.plotlist['U']
        x_CFD = CFDrhoObj.grid[0]
        y_CFD = CFDrhoObj.grid[1]

        nx = x_CFD.shape[0]
        ny = y_CFD.shape[0]

        dx = np.diff(x_CFD)[0]
        dy = np.diff(y_CFD)[1]

        Lx_CFD = x_CFD.max()+dx/2.
        Ly_CFD = y_CFD.max()+dy/2.

    #Keep the same loaded MD field for timestepratio CFD steps
    timestepratio = 10
    ntimestep = rhoObj.maxrec*timestepratio

    if plot:
        #Setup figure
        fig, ax = plt.subplots(1,1)
        plt.ion()
        plt.show()
        n = 0
        #Colormap
        tcmap=plt.cm.RdYlBu_r
        first_time = True


    #Start coupled run
    comm = MPI.COMM_WORLD
    CPL = CPL()

    MD_COMM = CPL.init(CPL.MD_REALM)
    MD_Cart_COMM = MD_COMM.Create_cart([2, 1, 1])
    CPL.setup_md(MD_Cart_COMM, 
                 xyzL=[476.22031559046, 476.22031559046, 9.5244063118092], 
                 xyz_orig=[0.0, 0.0, 0.0])


    MD_rank = MD_COMM.Get_rank()
    MD_coords = MD_Cart_COMM.Get_coords(MD_rank)

    recv_array, send_array = CPL.get_arrays(recv_size=3, send_size=4)
    uwall = 0.; vwall = 0.5
    olap_limits = np.zeros(6); portion = np.zeros(6)
    olap_limits = CPL.get_olap_limits()
    portion = CPL.my_proc_portion(olap_limits)
    dV = CPL.get("dx")*CPL.get("dy")*CPL.get("dz")

    mdrec = startrec #-1
    for timestep in range(1,ntimestep):

        if MD_rank == 0:
            print("CFD timestep = ", timestep, "MD record loaded=", mdrec, flush=True)

        if timestep%timestepratio == 0:
            mdrec = mdrec + 1

        recv_array, ierr = CPL.recv(recv_array)

        #Get density bottom boundary condition
        MDrho = rhoObj.read(startrec=mdrec, endrec=mdrec)
        MDrhoBC = np.mean(MDrho[:,yloc,:,:,:],(2,3))

        MDu = uObj.read(startrec=mdrec, endrec=mdrec)
        MDuBC = np.mean(MDu[:,yloc,:,:,:],(2))

        Npercell = MDrhoBC[portion[0]:portion[1]+1,portion[4]:portion[5]+1]*dV
        send_array[3,:,0,:] = Npercell

        uvw = MDuBC[portion[0]:portion[1]+1,portion[4]:portion[5]+1,:]
        send_array[0,:,0,:] = uvw[:,:,0]*Npercell
        send_array[1,:,0,:] = uvw[:,:,1]*Npercell
        send_array[2,:,0,:] = uvw[:,:,2]*Npercell

        CPL.send(send_array)

        if plot and MD_rank == 0:

            #PLOT MD density contour
            cm = ax.imshow(np.mean(MDrho[:,yloc:0:-1,:,:,:],(2,3,4)).T,
                           extent=[0.,Lx_MD,-yloc*dy_MD,0.],
                           cmap=tcmap, alpha=.5)

            #Plot CFD density contour
            if plot_CFD:
                CFDrho = CFDrhoObj.read(timestep, timestep)
                cm = ax.imshow(np.mean(CFDrho[:,::-1,:,:,:],(2,3,4)).T,
                           extent=[0.,Lx_CFD,0.,Ly_CFD],
                           cmap=tcmap)

            #Plot VMD DATA
            if plot_vmd:
                vmd = DCDObj[mdrec][:]
                vmd[:,0] = vmd[:,0]+.5*Lx_MD
                vmd[:,1] = vmd[:,1]+.5*Ly_MD - yloc*dy_MD
                mask = 0.>vmd[:,1]
                ax.plot(vmd[mask,0],vmd[mask,1],'ko',ms=0.1)

            #Plot CFD density, pressure or temperature
            if plot_grid:
                draw_grid(ax,int(nx/1.),int(ny/1.),nz=1, px=2,py=1,pz=1,
                          xmin=0.,ymin=0.,zmin=0.,xmax=Lx_CFD,ymax=Ly_CFD,zmax=1.,
                          lc='k',fc='k',label="CFD", XKCD_plots=False)

                draw_grid(ax, 1, 1, nz=1, px=2, py=1, pz=1,
                          xmin=0.,ymin=-yloc*dy_MD,zmin=0.,xmax=Lx_MD,ymax=0.,zmax=1.,
                          lc='k',fc='k',label="MD", XKCD_plots=False)
            if plot_quiver:
                if plot_CFD:
                    X, Y, UV = CFDUObj.quiver([0,1], startrec=timestep, endrec=timestep)
                    ax.quiver(X, Y, UV[:,:,0], UV[:,:,1])

            if first_time:
                plt.colorbar(cm)
                first_time=False

            ax.set_aspect('equal', 'box')
            #plt.axis("tight")
            plt.pause(0.01)
            if savefig:
                plt.savefig("cpl_bubble{0:04d}".format(n)+".png")
            plt.cla()
            n += 1

    #--------------------------#

