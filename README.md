# lowe-g4
Geant4 version of the low-energy electron scattering model.

# Installation
On UNIX systems:
Make sure you have CMake (3.3 or higher) installed
Install GDML and Xerces-C (see also the GEANT4 installation manual)
CLHEP is optional, GEANT4 has its own version included
Install GEANT4.10.02 patch 02 with GDML and Xerces-C enabled, a detailed description of this process can be found in the GEANT4 installation manual.
The GEANT4 DATA is also necessary, so be sure to install this during the GEANT4 installation.

Before installation, make sure to run the geant4.sh script to setup the right environment variables (this script is placed in the /bin/ in the install directory)
Install the low-energy extension by going into the CADPhysics directory
cd CADPhysics
and executing:
cmake .
make
make install

repeat this in the SEM directory:
cd ../SEM
cmake .
make
make install

the executable SEM_4.10.00 should now be placed inside the usual installation directory:
/usr/local/bin
or in the /bin/ directory in the specified installation directory

To run the example, make a symbolic link to the SEM_4.10.00 in the Examples/Transmission directory
cd ../Examples/Transmission
ln -s /usr/local/bin/SEM_4.10.00 SEM
now simply execute:
./run.sh
to run the transmission example
