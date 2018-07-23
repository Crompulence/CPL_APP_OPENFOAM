# start from cpl library base
FROM cpllibrary/cpl-library
MAINTAINER Edward Smith <edward.smith05@imperial.ac.uk>

#Number of processes to use in build
ENV FOAM_VERSION=3.0.1 \
    FOAM_INST_DIR="/OpenFOAM/"

ENV FOAM_SRC_DIR=$FOAM_INST_DIR/OpenFOAM-$FOAM_VERSION \
    APP_DIR=$FOAM_INST_DIR/CPL_APP_OPENFOAM-$FOAM_VERSION \
    WM_NCOMPPROCS=12
    

#Get OpenFOAM and APP
WORKDIR $FOAM_INST_DIR
RUN wget http://downloads.sourceforge.net/foam/OpenFOAM-$FOAM_VERSION.tgz && \
    tar -xvf OpenFOAM-$FOAM_VERSION.tgz && \
    rm -f ./OpenFOAM-$FOAM_VERSION.tgz && \
    wget http://downloads.sourceforge.net/foam/ThirdParty-$FOAM_VERSION.tgz && \
    tar -xvf ThirdParty-$FOAM_VERSION.tgz && \
    rm -f ./ThirdParty-$FOAM_VERSION.tgz && \
    git clone https://github.com/Crompulence/CPL_APP_OPENFOAM-$FOAM_VERSION.git $APP_DIR

#Get Prerequists
RUN apt-get update && apt-get install -y \
    bison \
    flex-old \
    libboost-system-dev \
    libboost-thread-dev \
    libncurses-dev  \
    libreadline-dev\
    libxt-dev \
    libz-dev \
 && rm -rf /var/lib/apt/lists/*

#We copy this pref file to build OpenFOAM with system MPICH instead of OpenMPI
RUN cp $APP_DIR/config/prefs_system_mpich.sh $FOAM_SRC_DIR/etc/prefs.sh

#Build from CPL APP file
RUN echo $FOAM_INST_DIR > $APP_DIR/CODE_INST_DIR

#Ideally we set all these as ENV but so many for OpenFOAM
WORKDIR $APP_DIR
RUN /bin/bash -c "source SOURCEME.sh && \
    cd $FOAM_INST_DIR/ThirdParty-$FOAM_VERSION && \
    ./Allwmake"

RUN /bin/bash -c "source SOURCEME.sh && \
    cd $FOAM_SRC_DIR && \
    mkdir -p platforms/linux64GccDPInt32OptSYSTEMMPI/src/Pstream/mpi && \
    mkdir -p platforms/linux64GccDPInt32OptSYSTEMMPI/src/parallel/decompose/ptscotchDecomp && \
    patch src/parallel/decompose/ptscotchDecomp/ptscotchDecomp.C $APP_DIR/config/ptscotchDecomp.patch && \
    ./Allwmake -j"

#We need to make Pstream, patch the OpenFOAM version and build solvers
RUN /bin/bash -c "source SOURCEME.sh && \
    make clean && \ 
    make pstream && \
    mv $FOAM_SRC_DIR/platforms/linux64GccDPInt32Opt/lib/mpi-system/libPstream.so $FOAM_SRC_DIR/platforms/linux64GccDPInt32Opt/lib/mpi-system/libPstream.so.orig && \
    cp lib/libPstream.so $FOAM_SRC_DIR/platforms/linux64GccDPInt32Opt/lib/mpi-system/ && \
    make"

 
