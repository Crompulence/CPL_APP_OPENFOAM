# Start from a fully built version of OpenFOAM and just rebuild the APP
# This is needed as OpenFOAM takes 8+ hours to build and dockerhub
# does not seem to allow caching
FROM cpllibrary/cplopenfoam
MAINTAINER Edward Smith <edward.smith05@imperial.ac.uk>

#We need to update the app to the latest code, patch OpenFOAM and build solvers
WORKDIR $APP_DIR
RUN git pull
RUN /bin/bash -c "source SOURCEME.sh && \
    make clean && \ 
    make pstream && \
    cp lib/libPstream.so $FOAM_SRC_DIR/platforms/linux64GccDPInt32Opt/lib/mpi-system/ && \
    make"

