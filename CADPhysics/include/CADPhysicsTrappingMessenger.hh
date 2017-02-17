// CADPhysicsTrappingMessenger.hh
//

#ifndef CADPhysicsTrappingMessenger_h
#define CADPhysicsTrappingMessenger_h 1

class CADPhysicsTrapping;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsTrappingMessenger: public G4UImessenger
{
public:
	CADPhysicsTrappingMessenger(CADPhysicsTrapping* mpga);
	~CADPhysicsTrappingMessenger();

public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
//	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsTrapping* target;

	G4UIdirectory*          trappingDir;

	// The only user command
	G4UIcmdWithADoubleAndUnit*	trap_C_Cmd;
	G4UIcmdWithADoubleAndUnit*	trap_g_Cmd;
};

#endif
