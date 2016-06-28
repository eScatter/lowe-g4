// CADPhysicsPhysicsList.hh
//
// This class is used to define which physical processes are used in the simulation, and their relative order is set. Also in some
// cases the desired behavior of the processes is specified.
// In the user's main file, an instance of this class is sent as an argument to the run manager's SetUserInitialization method.
//
// This particular physics list contains the 'standard' CADPhysics processes. It defines the following particles:
// - Geantino: a ficticious massless and chargeless particle without any interactions, used for debugging purposes;
// - Electron
// - Positron
// - Gamma
//
// Positrons and gammas are required for proper behavior of the LowEnergy processes, but gammas are also of interest in their
// own right.
// 
// The following processes are defined as acting on these particles:
// - For geantinos and positrons: only transportation
// - For electrons:
//		* Transportation; a special version that works together with the elastic scattering process;
//		* Elastic scattering for both solids and gases;
//		* Inelastic scattering for solids, following a dielectric function approach;
//		* A boundary crossing process, that takes the potential energy step at material boundaries into account and handles 
//		  transmission and reflection of electrons;
//		* Inelastic scattering (ionization) in gases: for this a modified version of the Geant4 low energy ionization process is used;
//		* Bremsstrahlung
//		* An artificial 'kill' process, that kills particles after they have travelled a given distance
// - For gammas:
//		* (Standard) transportation;
//		* Four Geant4 'low energy' processes for gammas. These model the photoelectric effect (photoionization), Compton scattering (scattering
//		  of gammas on electrons in the material), gamma conversion (into an electron/positron pair), and Rayleigh scattering (interaction
//		  with atoms).
//		  We have chosen to use the Penelope implementations of these processes.
//		  
// This physics list is considered suitable for most simulations of electron interaction with matter. Because the generation of gammas from 
// electrons and gamma processes are included, it can also be used for applications like EDX/X-ray microanalysis.
//
#ifndef CADPhysicsPhysicsList_h
#define CADPhysicsPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "CADPhysicsUnits.hh"

class CADPhysicsPhysicsList: public G4VUserPhysicsList
{
public:
	CADPhysicsPhysicsList();
	~CADPhysicsPhysicsList();

protected:
	// Construct particle(s) and physics
	// See G4VUserPhysicsList.hh:
	// ConstructParticle() is invoked directly by the RunManager, while
	// ConstructProcess() is invoked by the Construct method of G4VUserPhysicsList.
	void ConstructParticle();
	void ConstructProcess();

	void SetCuts();

};

#endif
