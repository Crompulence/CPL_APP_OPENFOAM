
# CPL APP for OpenFOAM [![Build Status](https://travis-ci.org/Crompulence/CPL_APP_OPENFOAM-3.0.1.svg?branch=master)](https://travis-ci.org/Crompulence/CPL_APP_OPENFOAM-3.0.1/) ![Build Status](https://img.shields.io/docker/cloud/build/cpllibrary/cpl-openfoam)

1 ) Pre-requisites for compilation
=================================


- A compiled version of CPL library (https://github.com/Crompulence/cpl-library) which has been added to the user's path using "source SOURCEME.sh"
- A copy of OpenFOAM
- The same GCC and MPICH version used to build CPL library

The following environnment variables:

    FOAM_CPL_SOCKET_LIBBIN = $FOAM_INST_DIR/cpl-socket/lib
    LD_LIBRARY_PATH        = $LD_LIBRARY_PATH:$FOAM_CPL_LIBBIN

**must be defined** in order for a) the compilation to work and b) the library
to be found by the ld linker. They are conveniently defined in the config 
file SOURCEME, which may be found in the level above this README file: 

    $  cd ../
    $  source SOURCEME.sh
    $  cd -

2 ) Install
===========

The install process assumes that you have OpenFOAM fully installed with all code source. The process then simply builds a top level solver, using all of the OpenFOAM shared libraries and CPL library's shared library, which is designed for coupled simulation. 

First, change directory to CPL_APP_OPENFOAM-DEV,

    cd CPL_APP_OPENFOAM-DEV

 and create a file called CODE_INST_DIR which will tell the APP which OPENFOAM to install against, 

    echo "/path/to/openfoam/directory/" > CODE_INST_DIR

Note that the folder in this directory is expected to be called OpenFOAM-3.0.1. 
This is the currently supported OpenFOAM version.

Now, to build the various solvers, it may be as simple as,

    $  make

N.B.: warnings from the included MPI headers may be ignored. This makes a range of OpenFOAM solvers, including

   - CPLIcoFOAM - A basic fluid solver for MD-CFD coupling which partially overlaps an MD simulation
   - CPLSediFOAM - SediFOAM (https://github.com/xiaoh/sediFoam) adapted to use the CPL library philosophy of minimal linked library
   - CPLCFDDEMFOAM - CFDDEM (https://www.cfdem.com/) adapted to use the CPL library philosophy of minimal linked library

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
 - test - a range of test cases run automatically on Travis CI
 - config - scripts to specify version of MPI to build OpenFOAM with and patches for scotch.

New folders created by building process

 - lib - dynamic-link library binaries are created in a new folder `./lib/`. 
 - bin - Executable solver applications are created in the a folder `./bin`.








