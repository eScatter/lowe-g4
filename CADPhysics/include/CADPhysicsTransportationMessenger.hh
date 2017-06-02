// CADPhysicsTransportationMessenger.hh
//

#ifndef CADPhysicsTransportationMessenger_h
#define CADPhysicsTransportationMessenger_h 1

class CADPhysicsTransportation;
class G4UIdirectory;
class G4UIcmdWithABool;

#include "G4UImessenger.hh"
#include "globals.hh"

class CADPhysicsTransportationMessenger: public G4UImessenger
{
public:
	CADPhysicsTransportationMessenger(CADPhysicsTransportation* mpga);
	~CADPhysicsTransportationMessenger();

public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
//	G4String GetCurrentValue(G4UIcommand * command);

private:
	CADPhysicsTransportation* target;

	G4UIdirectory*          transportationDir;

	// The only user command
  G4UIcmdWithABool*           outputCmd;
};

#endif
