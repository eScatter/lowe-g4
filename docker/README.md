Docker image for Runnig Lowe<sup>-</sup>G4
==========================================

We tested the current code with Geant version 4.10.02.p03. To build the Docker image, download precisely this version of Geant,

        wget http://geant4.web.cern.ch/geant4/support/source/geant4.10.02.p02.tar.gz

or follow [this link](http://geant4.web.cern.ch/geant4/support/source/geant4.10.02.p02.tar.gz), and place the file here manually. Then build the image by running,

        docker build -t geant4 .

and get some coffee. We have a small example in `Examples/Transmission`. This should run without errors.

