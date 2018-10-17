import pytest
import os
import sys
import shutil
import numpy as np
import subprocess as sp
import itertools
from OpenFOAM_vs_analytical import check_OpenFOAM_vs_Analytical

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# Import symwraplib
sys.path.insert(0, "./SimWrapPy/")
try:
    import simwraplib as swl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/SimWrapPy.git ./SimWrapPy"
    downloadout = sp.check_output(cmd, shell=True)
    print(downloadout)
    sys.path.insert(0, "./SimWrapPy")
    import simwraplib as swl

#Define test directory based on script file
TEST_DIR = os.path.dirname(os.path.realpath(__file__))

#Get all permutations
nx = 6; nz = 6
px = 1; pz = 1
params = []
for i in list(itertools.combinations_with_replacement(range(1,3), 4)):
    params.append([nx*i[2], nz*i[3], px*i[0], pz*i[1]])


@pytest.mark.parametrize("nx, nz, px, pz", params)
def test_newtest(nx, nz, px, pz):

    # Inputs that are the same for every thread
    basedir = TEST_DIR
    srcdir = None
    executable = "CPLIcoFoam"
    inputfile = "/openfoam_ico"
    rundir = TEST_DIR + "/run" + str(nx) + "_" + str(nz) + "_" + str(px) + "_" + str(pz)

    #Clean previous result, generate grid and decompose for parallel run
    with cd (TEST_DIR+"/"+inputfile):
        sp.check_output("python clean.py -f", shell=True)
        sp.check_output("blockMesh", shell=True)
        sp.check_output("decomposePar", shell=True)

    #Setup Changes
    with cd(TEST_DIR):

        changes={'cell':[nx,8,nz], 'process':[px,1,pz]}


        #Setup a LAMMPS run object
        of = swl.OpenFOAMRun(None, basedir, rundir,
                             executable, inputfile,
                             inputchanges=changes)

        #Setup a mock script
        mockscript = "./python_dummy/MD_sendrecv.py"
        mock = swl.ScriptRun(rundir, mockscript, 
                             inputchanges={"npxyz =": "["+str(px)+", 1, "+str(pz)+"]"})
        mock.get_nprocs = lambda : px*pz

        #Setup a coupled run
        topo = [1, nx, 1, 4, 1, nz]
        run = swl.CPLRun(None, basedir, rundir, [mock, of],
                         inputfile="/cpl/COUPLER.in", 
                         inputchanges={"OVERLAP_EXTENTS" : topo, 
                                       "CONSTRAINT_INFO" : [3, 0] + topo, 
                                       "BOUNDARY_EXTENTS": topo})
        #Run the case
        run.setup()
        run.execute(blocking=True, print_output=True, extra_cmds="-Mp")

        #Check results are correct
        check_OpenFOAM_vs_Analytical(rundir)

