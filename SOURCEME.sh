#!/bin/bash
#~~~
#    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
#     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
#      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
#       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
#        _\/\\\_____________\/\\\/////////____\/\\\_____________
#         _\//\\\____________\/\\\_____________\/\\\_____________
#          __\///\\\__________\/\\\_____________\/\\\_____________
#           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
#            _______\/////////__\///______________\///////////////__
#~~~
#

# Environment variable for install directory
export FOAM_CPL_VERSION=3.0.1
export CWD=`pwd`
export FOAM_INST_DIR=$(cat CODE_INST_DIR)
# Source the other environment variables
foamDotFile=$FOAM_INST_DIR/OpenFOAM-$FOAM_CPL_VERSION/etc/bashrc
echo $foamDotFile
if [ -f $foamDotFile ]; then
    source $foamDotFile
else
    echo "ERROR:"
    echo "   Configuration file 'OpenFOAM-$FOAM_CPL_VERSION/etc/bashrc' not found."
    return 1
fi

echo ""
echo "FOAM_MPI environment variable is now: " $FOAM_MPI
export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$FOAM_MPI

# CPL APP environment variables
echo ""
echo "New environment variables: "
echo ""

export FOAM_CPL_APP=$CWD
echo "    FOAM_CPL_APP = " $FOAM_CPL_APP

export FOAM_CPL_APP_SRC=$FOAM_CPL_APP/src
echo "    FOAM_CPL_APP_SRC = " $FOAM_CPL_APP_SRC

export FOAM_CPL_APP_LIBBIN=$FOAM_CPL_APP/lib
echo "    FOAM_CPL_APP_LIBBIN = " $FOAM_CPL_APP_LIBBIN

export FOAM_CPL_APP_BIN=$FOAM_CPL_APP/bin
echo "    FOAM_CPL_APP_BIN = " $FOAM_CPL_APP_BIN

echo ""

# Update paths
export PATH=$PATH:$FOAM_CPL_APP_BIN
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FOAM_CPL_APP_LIBBIN
echo "PATH updated to:"
echo "   $PATH:"$FOAM_CPL_APP_BIN
echo
#echo "LD_LIBRARY_PATH updated to: "
#echo "   $LD_LIBRARY_PATH:"$FOAM_CPL_APP_LIBBIN
