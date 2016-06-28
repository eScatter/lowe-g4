// CADPhysicsIonBoundaryMessenger.hh
//

#ifndef CADPhysicsIonBoundaryMessenger_h
#define CADPhysicsIonBoundaryMessenger_h 1

class CADPhysicsIonBoundary;
class G4UIdirectory;
class G4UIcmdWithADouble;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsIonBoundaryMessenger: public G4UImessenger
{
public:
	CADPhysicsIonBoundaryMessenger(CADPhysicsIonBoundary* mpga);
	~CADPhysicsIonBoundaryMessenger();

	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsIonBoundary* target;

	G4UIdirectory*          IBDir;

	// The only user command
	G4UIcmdWithADouble*		gammaCmd;
};

#endif
