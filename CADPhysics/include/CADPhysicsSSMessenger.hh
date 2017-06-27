// CADPhysicsSSMessenger.hh
//

#ifndef CADPhysicsSSMessenger_h
#define CADPhysicsSSMessenger_h 1

class CADPhysicsSingleScattering;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsSSMessenger: public G4UImessenger
{
public:
	CADPhysicsSSMessenger(CADPhysicsSingleScattering* mpga);
	~CADPhysicsSSMessenger();

	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsSingleScattering* target;

	G4UIdirectory*          SSDir;

	// User commands
	G4UIcmdWithABool*		multistepCmd;
	G4UIcmdWithAnInteger*	verboseCmd;

};

#endif
