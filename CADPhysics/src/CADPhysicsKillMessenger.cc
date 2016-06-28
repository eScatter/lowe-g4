// CADPhysicsKillMessenger.cc
//

#include "CADPhysicsKillMessenger.hh"
#include "CADPhysicsKill.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

CADPhysicsKillMessenger::CADPhysicsKillMessenger(CADPhysicsKill * mpga)
:target(mpga)
{
	killDir = new G4UIdirectory("/process/kill/");
	killDir->SetGuidance("Controls for the CADPhysicsKill process.");

	lengthCmd = new G4UIcmdWithADoubleAndUnit("/process/kill/length",this);
	lengthCmd->SetGuidance("Set the total distance electrons may travel before being killed.");
	lengthCmd->SetParameterName("length",true);
	lengthCmd->SetDefaultValue(1);
	lengthCmd->SetDefaultUnit("mm");
	lengthCmd->SetUnitCandidates("nanometer micron mm cm m km");
	lengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

CADPhysicsKillMessenger::~CADPhysicsKillMessenger()
{
	delete lengthCmd;
	delete killDir;
}

void CADPhysicsKillMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
	if( command==lengthCmd )
	{ target->SetLength(lengthCmd->GetNewDoubleValue(newValue)); }
}

G4String CADPhysicsKillMessenger::GetCurrentValue(G4UIcommand * command)
{
	// Print the 'kill length' by invoking '?process/kill/length' from the command line.
	G4String cv;
	if( command==lengthCmd )
	{ cv = lengthCmd->ConvertToString(target->GetLength()); }
	return cv;
}
