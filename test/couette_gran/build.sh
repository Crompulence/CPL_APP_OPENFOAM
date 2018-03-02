#! /bin/bash
MYPWD=${PWD}
ROOTDIR=../../../
CPL_DIR=/home/es205/codes/cpl/cpl-library/
LAMMPS_DIR=${ROOTDIR}/LAMMPS-dev_coupled/CPL_APP_LAMMPS-DEV/
OPENFOAM_DIR=${ROOTDIR}/OpenFOAM-3.0.1_coupled/CPL_APP_OPENFOAM-3.0.1/

if [ "$1" -eq "0" ]
then
    cd ${CPL_DIR}
    make
    cd ${MYPWD}

    rm -f ./lmp_cpl
    cd ${LAMMPS_DIR}
    make 
    cp ./bin/lmp_cpl ${MYPWD}
    cd ${MYPWD}
    cd ${MYPWD}

    rm -f openfoam/CPLSediFOAM
    cd ${ROOTDIR}/OpenFOAM-3.0.1_coupled/CPL_APP_OPENFOAM-3.0.1/
    make sedifoam 2>&1 | grep -i 'error' 
    cp ./bin/CPLSediFOAM ${MYPWD}/openfoam/
    cd ${MYPWD}
fi

#Build CPL library
if [ "$1" -eq "1" ]
then
    cd ${CPL_DIR}
    make
    cd ${MYPWD}
fi

#Build LAMMPS
if [ "$1" -eq "2" ]
then
    rm -f ./lmp_cpl
    cd ${LAMMPS_DIR}
    make 
    cp ./bin/lmp_cpl ${MYPWD}
    cd ${MYPWD}
    cd ${MYPWD}
fi

#Build OpenFOAM
if [ "$1" -eq "3" ]
then
    rm -f openfoam/CPLSediFOAM
    cd ${OPENFOAM_DIR}
    make sedifoam 2>&1 | grep -i 'error' 
    cp ./bin/CPLSediFOAM ${MYPWD}/
    cd ${MYPWD}
fi


