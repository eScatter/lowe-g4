// CADPhysicsDIMessenger.hh
//

#ifndef CADPhysicsDIMessenger_h
#define CADPhysicsDIMessenger_h 1

class CADPhysicsDI;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsDIMessenger: public G4UImessenger
{
public:
	CADPhysicsDIMessenger(CADPhysicsDI* mpga);
	~CADPhysicsDIMessenger();

	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsDI* target;

	G4UIdirectory* DIDir;

	// User commands
	G4UIcmdWithABool*			secondariesCmd;
	G4UIcmdWithABool*			xrayCmd;
	G4UIcmdWithABool*			augerCmd;
	G4UIcmdWithABool*			rangecutCmd;
	G4UIcmdWithAnInteger*		verboseCmd;
	G4UIcmdWithADoubleAndUnit*	energycutCmd;
	G4UIcommand*				resetcounterCmd;
  G4UIcmdWithABool*           outputCmd;
};

#endif
