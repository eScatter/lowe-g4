// CADPhysicsLowEnIoniMessenger.hh
//

#ifndef CADPhysicsLowEnIoniMessenger_h
#define CADPhysicsLowEnIoniMessenger_h 1

class CADPhysicsLowEnergyIonisation;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsLowEnIoniMessenger: public G4UImessenger
{
public:
	CADPhysicsLowEnIoniMessenger(CADPhysicsLowEnergyIonisation* mpga);
	~CADPhysicsLowEnIoniMessenger();

public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsLowEnergyIonisation* target;

	G4UIdirectory*				LEIDir;

	// User commands
	G4UIcmdWithABool*			usecutCmd;
	G4UIcmdWithADoubleAndUnit*	ecutCmd;
	G4UIcmdWithABool*			generateionsCmd;
};

#endif
