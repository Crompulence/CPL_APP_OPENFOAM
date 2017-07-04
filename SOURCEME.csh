#!/bin/csh
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
setenv FOAM_CPL_VERSION 3.0.1
setenv CWD=`pwd`
setenv FOAM_INST_DIR=$(cat CODE_INST_DIR)
# Source the other environment variables
set foamDotFile=$FOAM_INST_DIR/OpenFOAM-$FOAM_CPL_VERSION/etc/cshrc
if( -f $foamDotFile ) then
    source $foamDotFile
else
    echo "ERROR:"
    echo "   Configuration file 'OpenFOAM-$FOAM_CPL_VERSION/etc/cshrc' not found."
    return 1
endif

echo ""
echo "FOAM_MPI environment variable is now: " $FOAM_MPI
setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/$FOAM_MPI

# CPL APP environment variables
echo ""
echo "New environment variables: "
echo ""

setenv FOAM_CPL_APP $CWD
echo "    FOAM_CPL_APP = " $FOAM_CPL_APP

setenv FOAM_CPL_APP_SRC $FOAM_CPL_APP/src
echo "    FOAM_CPL_APP_SRC = " $FOAM_CPL_APP_SRC

setenv FOAM_CPL_APP_LIBBIN $FOAM_CPL_APP/lib
echo "    FOAM_CPL_APP_LIBBIN = " $FOAM_CPL_APP_LIBBIN

setenv FOAM_CPL_APP_BIN $FOAM_CPL_APP/bin
echo "    FOAM_CPL_APP_BIN = " $FOAM_CPL_APP_BIN

echo ""

# Update paths
setenv PATH $PATH\:$FOAM_CPL_APP_BIN
setenv LD_LIBRARY_PATH $FOAM_CPL_APP_LIBBIN:$LD_LIBRARY_PATH
echo "PATH updated to:"
echo "   $PATH\:$FOAM_CPL_APP_BIN"
echo
echo "LD_LIBRARY_PATH updated to:"
echo "   $FOAM_CPL_APP_LIBBIN:$LD_LIBRARY_PATH"
