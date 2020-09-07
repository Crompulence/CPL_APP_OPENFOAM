# Start from a fully built version of OpenFOAM and just rebuild CPL library
# and the interface to this in the form of CPL_APP_OPENFOAM
# This is needed as OpenFOAM takes 8+ hours to build and dockerhub
# does not seem to allow caching

FROM cpllibrary/mpich-openfoam
MAINTAINER Edward Smith <edward.smith05@imperial.ac.uk>

# Check if CPL library needs to be updated, this is required as cplopenfoam 
# will contain an old version of CPL library and this can easily be updated
WORKDIR /cpl-library
RUN git fetch --all && \
    git reset --hard origin/master && \
    make PLATFORM=gcc clean && \ 
    make PLATFORM=gcc

#We need to update the app to the latest code, patch OpenFOAM and build solvers
WORKDIR $APP_DIR
RUN git fetch --all && \
    git reset --hard origin/master
RUN /bin/bash -c "source SOURCEME.sh && \
    make clean && \ 
    make pstream && \
    cp lib/libPstream.so $FOAM_SRC_DIR/platforms/linux64GccDPInt32Opt/lib/mpi-system/ && \
    make"

