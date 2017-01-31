FROM debian:8.7
MAINTAINER Johan Hidding <j.hidding@esciencecenter.nl>

RUN echo "deb http://ftp.debian.org/debian jessie-backports main" \
        >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install --no-install-recommends -y \
        build-essential libxerces-c-dev
RUN apt-get -t jessie-backports install --no-install-recommends -y \
        cmake
RUN apt-get install -y \
        libexpat1-dev git vim screen

WORKDIR /usr/src
COPY geant4.10.02.p03.tar.gz .
RUN tar xf geant4.10.02.p03.tar.gz
RUN        mkdir build \
        && cd build \
        && cmake ../geant4.10.02.p03 \
                -DGEANT4_INSTALL_DATA=ON \
                -DGEANT4_BUILD_MULTITHREADED=ON \
                -DGEANT4_USE_GDML=ON \
        && make -j4 \
        && make install

# Bug in Geant4 install script
RUN cp geant4.10.02.p03/source/visualization/VRML/include/G4FRClient.hh /usr/local/include/Geant4/

RUN /bin/bash -c "source /usr/local/bin/geant4.sh"

RUN git clone http://github.com/eScatter/lowe-g4.git
RUN        mkdir build-CAD \
        && cd build-CAD \
        && cmake ../lowe-g4/CADPhysics \
        && make && make install

RUN rm lowe-g4/SEM/src/G4PolyconeSide.cc
RUN        mkdir build-SEM \
        && cd build-SEM \
        && cmake ../lowe-g4/SEM \
        && make && make install

