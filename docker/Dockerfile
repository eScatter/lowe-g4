FROM debian:8.7
MAINTAINER Johan Hidding <j.hidding@esciencecenter.nl>

RUN echo "deb http://ftp.debian.org/debian jessie-backports main" \
        >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install --no-install-recommends -y \
        build-essential libxerces-c-dev wget
RUN apt-get -t jessie-backports install --no-install-recommends -y \
        cmake
RUN apt-get install -y \
        libexpat1-dev git vim screen

WORKDIR /usr/src
COPY geant4.10.02.p02.tar.gz .
RUN tar xf geant4.10.02.p02.tar.gz
RUN        mkdir build \
        && cd build \
        && cmake ../geant4.10.02.p02 \
                -DGEANT4_INSTALL_DATA=ON \
                -DGEANT4_BUILD_MULTITHREADED=ON \
                -DGEANT4_USE_GDML=ON \
        && make -j4 \
        && make install

# Bug in Geant4 install script
RUN cp geant4.10.02.p02/source/visualization/VRML/include/G4FRClient.hh /usr/local/include/Geant4/

RUN /bin/bash -c "source /usr/local/bin/geant4.sh"

RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz \
 && tar xf hdf5-1.10.1.tar.gz \
 && cd hdf5-1.10.1 \
 && ./configure --enable-cxx --prefix=/usr/local \
 && make && make install

RUN git clone -b hdf5 https://github.com/eScatter/lowe-g4.git \
 && cd lowe-g4 && git submodule init && git submodule update

RUN        mkdir build-CAD \
        && cd build-CAD \
        && cmake ../lowe-g4/CADPhysics \
        && make && make install

RUN        mkdir build-SEM \
        && cd build-SEM \
        && cmake ../lowe-g4/SEM \
        && make && make install
        
RUN        ln -s /usr/local/bin/SEM_4.10.00 /usr/local/bin/SEM
RUN        echo "source /usr/local/bin/geant4.sh" >> /etc/bash.bashrc
