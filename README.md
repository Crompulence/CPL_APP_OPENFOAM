
# CPL APP for OpenFOAM [![CI](https://github.com/Crompulence/CPL_APP_OPENFOAM/actions/workflows/main.yml/badge.svg)](https://github.com/Crompulence/CPL_APP_OPENFOAM/actions/workflows/main.yml)

1 ) Pre-requisites for compilation
=================================


- A compiled version of CPL library (https://github.com/Crompulence/cpl-library) which has been added to the user's path using "source SOURCEME.sh"
- A copy of OpenFOAM
- The same GCC and MPICH version used to build CPL library

The following environment variables:

    FOAM_CPL_SOCKET_LIBBIN = $FOAM_INST_DIR/cpl-socket/lib
    LD_LIBRARY_PATH        = $FOAM_CPL_LIBBIN:$LD_LIBRARY_PATH

**must be defined** in order for a) the compilation to work and b) the library
to be found by the ld linker. They are conveniently defined in the config 
file SOURCEME.sh, located in cpl-library, which can be used as follows: 

    $  source SOURCEME.sh

2 ) Install
===========

The install process assumes that you have OpenFOAM fully installed with all source code. The process then simply builds a top level solver, using all of the OpenFOAM shared libraries and CPL library's shared library, which is designed for coupled simulation. 

First, change directory to CPL_APP_OPENFOAM,

    cd CPL_APP_OPENFOAM

 and create a file called CODE_INST_DIR which will tell the APP which OPENFOAM to install against.  (OpenFOAM typically makes this path available via the environment variable $FOAM_INSTALL_DIR) 

    echo "/path/to/openfoam/directory/" > CODE_INST_DIR

Note that the folder in this directory is expected to be called OpenFOAM-v2106. Currently supported OpenFOAM versions include v2106, v2112 and previously 3.0.1. In addition, a version using [Foam-extend](https://github.com/FoamScience/CPL_FoamExtend) has been developed.

Next, source the SOURCEME.sh file. NB this is the SOURCEME.sh file which sits inside the CPL_APP_OPENFOAM directory and is different to the SOURCEME.sh file which sits inside the cpl-library directory.

    $  source SOURCEME.sh

Now, to build the various solvers, it may be as simple as,

    $  make

N.B.: warnings from the included MPI headers may be ignored. This makes a range of OpenFOAM solvers, including

   - CPLTestFoam - A version of an OpenFOAM solver with the same contents as the [CPL library minimal send and recieve mocks](https://github.com/Crompulence/cpl-library/tree/master/examples/minimal_send_recv_mocks), this tests to see if both CPL library and OpenFOAM is built correctly and exchnages information using the code in examples/CPLTestFoam. 
   - CPLTestSocketFoam - A minimal test of the functionality of the CPLSocket code which is used to exchange information and set boundary conditions in OpenFOAM, run using the examples/CPLTestSocketFoam with debugging information shown.
   - CPLIcoFOAM - A basic fluid solver for MD-CFD coupling which partially overlaps an MD simulation, same basic code as CPLTestSocketFoam but without the degging information.
   - CPLinterCondensatingEvaporatingFoam - Used for multi-phase simulations and runs with the example in examples/interCondensatingEvaporatingFoam.
   - CPLSediFOAM - SediFOAM (https://github.com/xiaoh/sediFoam) adapted to use the CPL library philosophy of minimal linked library.
   - CPLCFDDEMFOAM - CFDDEM (https://www.cfdem.com/) adapted to use the CPL library philosophy of minimal linked library.

3 ) License
==========

CPL_APP_OpenFOAM is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.


4 ) Directory Structure
=========================

This application repository is structured as follows:

 - src - source files which include the following
   - CPLPstream: OpenFOAM's parallel-computing communications supports a number of different
paradigms, including MPI. In order to support all of them, it wraps all 
communication functions in a library called "Pstream". CPL Library is based
on MPI, but requires the MPI_COMM_WORLD communicator to be split into 
two "realm" communicators - one for MD, and the other for CFD. All 
subsequent operations that would have been on MPI_COMM_WORLD in the MD
and CFD coupled programs must now be called on the realm communicator
CPL_REALM_COMM. CPLPstream replicates OpenFOAM's MPI Pstream, where CPL_REALM_COMM has
been substituted in place of MPI_COMM_WORLD.
   - CPLSocket: The interface to CPL Library is implemented in a single class named
CPLSocket. This directory holds its source code. 
   - solvers: OpenFOAM is implemented as a set of libraries, and it is up to the 
user to develop their own main-level applications that use these 
libraries. This directory contains the source code for a coupled
incompressible solver CPLIcoFoam that is based on OpenFOAM's icoFoam.

 - examples - some input examples for various cases
 - test - a range of test cases run automatically on GitHub Actions (previously Travis CI)
 - config - scripts to specify version of MPI to build OpenFOAM (in particular mpich) as well as other patches.

New folders created by building process

 - lib - dynamic-link library binaries are created in a new folder `./lib/`. 
 - bin - Executable solver applications are created in the a folder `./bin`.



5 ) To patch a Pstream
=========================

Why do you need to do this? If you want to run two codes which share a single MPI_COMM_WORLD,

    mpiexec -n 4 CPLIcoFoam : -n 16 ./MD

which we will call the "shared" paradigm of coupling. This is as opposed to the distinct paradigm where both codes have their own distinct `MPI_COMM_WORLD` and are started individually then join together. The distinct paradigm relies on the ((not always functional) `MPI_Open_port` and `MPI_Comm_accept` linking to create an intercommunicator between the `MPI_COMM_WORLD` intracommunicators of both codes. The sharing of `MPI_COMM_WORLD` in shared coupling means that any use of `MPI_COMM_WORLD` in any MPI communications will now cause errors or deadlock in the coupled code, so these have to be replaced with a local comm (i.e. a CFD_WORLD_COMM). Patching OpenFOAMS's Pstream, the location where all MPI communication is contained, is the easiest way to do that for OpenFOAM. The steps are as follows (given in general terms to account for future OpenFOAM changes but specifically done up to v2112).

The Pstream which is used can be replaced for all codes using the `LD_LIBRARY_PATH` environment variable, such that the location of the Pstream library, namely libPstream.so, changes from a version with a path like

    /home/USERNAME/codes/CFD/OpenFOAM/openfoam-OpenFOAM-v2112.220610/platforms/linux64GccDPInt32Opt/lib/Pstream.so
    
to something like

    /home/USERNAME/CPL_APP_OPENFOAM/lib/Pstream.so
    
We will now describe the steps to patch Pstream where we check for correctness after each step.

The first step in patching requires recursively copying all of the Pstream/mpi file from the OpenFOAM code diretory. e.g. in my case

    /home/USERNAME/codes/CFD/openfoam-OpenFOAM-v2112.220610/src/Pstream/mpi/

to a folder called `CPLPstream` in the `CPL_APP_OPENFOAM/src` directory, replacing whatever files are there.

You can try building this with wmake, which should add PStream to `CPL_APP_OPENFOAM/lib` and
by changing `LD_LIBRARY_PATH` you should see included `CPL_APP_OPENFOAM/lib/Pstream.so`
replacing the default `Pstream.so` in your OpenFOAM executables.  For example, the following command displays that libraries are included in the icoFoam executable provided by your OpenFOAM, where you'll see your version of Pstream.so is now included.

    ldd /home/USERNAME/codes/CFD/OpenFOAM/openfoam-OpenFOAM-v2112.220610/platforms/linux64GccDPInt32Opt/bin/icoFoam

Now, you can edit this copy of Pstream in `CPL_APP_OPENFOAM/src/CPLPstream` as needed to ensure the shared CPL paradigm works. The following line must be added as the last decleration line near the bottom of the file `PstreamGlobals.H`:

    extern MPI_Comm CPLRealmComm;

and then setting CPLRealmComm to be `MPI_COMM_WORLD` for the default case which is nothing to do with coupling.  The file `PstreamGlobals.C` must now include the following line as the last decleration near the the bottom of the file:

    MPI_Comm Foam::PstreamGlobals::CPLRealmComm = MPI_COMM_WORLD; 

At this stage, you can rebuild Pstream again and functionality should be identical (we have just add a new variable so far).

Next, find and replace every instance of `MPI_COMM_WORLD` in `UPStream.C` with `Foam::PstreamGlobals::CPLRealmComm`. You can rebuild Pstream and test again as nothing has changed, as the new `Foam::PstreamGlobals::CPLRealmComm` is still `MPI_COMM_WORLD`.

Finally, in order to allow a "shared" MPI run, we need to define the `Foam::PstreamGlobals::CPLRealmComm` to be the value returned by `CPL.init`. In an example code this looks like:


    #include "PstreamGlobals.H"
    #include "mpi.h"
    #include "cpl.h"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    int main(int argc, char *argv[])
    {

        //Define variables
        int CFD_realm = 1;
        MPI_Comm CFD_COMM;

        int flag = 0;
        int ierr = MPI_Initialized(&flag);
        if (flag == 0)
		    MPI_Init(&argc, &argv);

        //Initialise CPL library
        CPL::init(CFD_realm, CFD_COMM); 
	    Foam::PstreamGlobals::CPLRealmComm = CFD_COMM;

A few notes, the PstreamGlobals.H must be included so the variable introduced above can be set (and replace `MPI_COMM_WORLD`). Also, this should be called as early as possible in the main function of an OpenFOAM solver. The function which starts MPI `Foam::"UPstream::init"` and uses `Foam::PstreamGlobals::CPLRealmComm` is called by one of these three include statements `#include "setRootCase.H"`, `#include "createTime.H"` or `#include "createMesh.H"`, so therefore the two lines 

    CPL::init(CFD_realm, CFD_COMM); 
    Foam::PstreamGlobals::CPLRealmComm = CFD_COMM;

must come immediately after `MPI_init`, and the `include PstreamGlobals.H` must come before any of these three include statements.








