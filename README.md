# lowe-g4
Geant4 version of the low-energy electron scattering model.

# Installation
On UNIX systems:
* Make sure you have CMake (3.3 or higher) installed
* Install GDML and Xerces-C (see also the GEANT4 installation manual)
* CLHEP is optional, GEANT4 has its own version included
* Install GEANT4.10.02 patch 02 with GDML and Xerces-C enabled, a detailed description of this process can be found in the GEANT4 installation manual.
* The GEANT4 DATA is also necessary, so be sure to install this during the GEANT4 installation.

Before installation, make sure to run the `geant4.sh` script to setup the right environment variables (this script is placed in the /bin/ in the install directory)

Install the low-energy extension by going into the CADPhysics directory `cd CADPhysics` and executing:

    cmake .
    make
    make install

repeat this in the SEM directory:

    cd ../SEM
    cmake .
    make
    make install

the executable SEM_4.10.00 should now be placed inside the usual installation directory `/usr/local/bin`
or in the `/bin/` directory in the specified installation directory.

To run the example, make a symbolic link to the SEM_4.10.00 in the `Examples/Transmission` directory

    cd ../Examples/Transmission
    ln -s /usr/local/bin/SEM_4.10.00 SEM

# Running a simulation
After the installation, simply execute the following command in the `Examples/Transmission`:
 
    ./run.sh

If the install directory was set to a different directory than the default directory, the environment variable `CADPHYSICS_BASE` has to be set to the install path `INSTALL_DIRECTORY`. In this case the example is started by:

    CADPHYSICS_BASE=INSTALL_DIRECTORY ./run.sh

This will start a simulation for alumina and a simulation for silicon as can be seen when `run.sh` is opened in a text editor.

# Setting up the simulation parameters
What exactly is simulated can be set in the files present in `Examples/Transmission`. The sample material can be set in `run.sh` as mentioned before.
In `YieldCurve.mac` the following can be set:

The amount of primary particles per energy. In the example this is set to 500 in order to have a quick run:

    /control/alias NPE 500

The detectors can be set directional or not, if they are directional, the direction has to be set.

    /detectors/PlaneR/directional true
    /detectors/PlaneR/direction 0.0 0.0 1.0

In the example, the sample itself is also set as a detector, in this case multiple hits is enabled, so that the particles can be tracked inside the sample:

    /detectors/Sample/multipleHits true

The primary gun is placed in the geometry:

    /gps/pos/centre 0.0 0.0 0.0001 mm

The material of the sample is set to the material given in `run.sh`:

    /detectors/setmaterial LogSample {material}
And at last the energies for which the yield is simulated are set. The energy has to be given in eV:

    /control/foreach mac/energy.mac energy 1250
Note that in the expample only one energy is given, if multiple energies are desired e.g. 500, 1000 and 1500 eV:

    /control/foreach mac/energy.mac energy 500 1000 1500

The simulation geometry is defined in `BulkSample.gdml` and `YieldCurve.gdml`.
In `BulkSample.gdml` first the solids are defined, then the logical volumes are defined using the already defined solids. Then the logical volumes are placed in the `samplecontainer` and the physical volumes are defined.
The geometry in the example consists of a 500 nm thick membrane as sample with two detectors. One detector is placed above the sample, the reflection detector and one detector is placed below the sample, the transmission detector. Note that the top of the sample is at 0.0 on the z-axis, which means the primary gun is positioned above the sample in `YieldCurve.mac`.
In `YieldCurve.gdml` the `samplecontainer` is placed inside the `World` volume and the list of detectors is defined.

In the `mac` directory, two more mac files are present, `Pointsource.mac` and `energy.mac`.
In `Pointsource.mac` the particle source is setup.
`energy.mac` sets up the energy of the gun, not that here the energies from `YieldCurve.mac` are used.
It is also possible to set an energy window of the detectors (commented out in the example)
The output files are also se in `energy.mac`. In the example, the yields of the reflection and transmission detector are set as output as well as the hits for both detectors and the sample.
The hits of the sample can be used to reconstruct all electron tracks. If this is not necessary, the line

    /detectors/Sample/outputHits {directory}/EmissionSample_{material}_{energy}.dat
can be commented out. In this case, the sample can also be removed from the detectorlist in `YieldCurve.mac`.
