// CADPhysicsKillMessenger.hh
//

#ifndef CADPhysicsKillMessenger_h
#define CADPhysicsKillMessenger_h 1

class CADPhysicsKill;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsKillMessenger: public G4UImessenger
{
public:
	CADPhysicsKillMessenger(CADPhysicsKill* mpga);
	~CADPhysicsKillMessenger();

public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsKill* target;

	G4UIdirectory*          killDir;

	// The only user command
	G4UIcmdWithADoubleAndUnit*	lengthCmd;
};

#endif


