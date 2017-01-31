// CADPhysicsTrappingMessenger.cc
//

#include "CADPhysicsTrappingMessenger.hh"
#include "CADPhysicsTrapping.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

CADPhysicsTrappingMessenger::CADPhysicsTrappingMessenger(CADPhysicsTrapping * mpga)
:target(mpga)
{
	trappingDir = new G4UIdirectory("/process/trapping/");
	trappingDir->SetGuidance("Controls for the CADPhysicsTrapping process.");

	//lengthCmd = new G4UIcmdWithADoubleAndUnit("/process/trapping/length",this);
	//lengthCmd->SetGuidance("Set the total distance electrons may travel before being trapped.");
	//lengthCmd->SetParameterName("length",true);
	//lengthCmd->SetDefaultValue(1);
	//lengthCmd->SetDefaultUnit("mm");
	//lengthCmd->SetUnitCandidates("nanometer micron mm cm m km");
	//lengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

CADPhysicsTrappingMessenger::~CADPhysicsTrappingMessenger()
{
	//delete lengthCmd;
	delete trappingDir;
}

//void CADPhysicsTrappingMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
//{
//	if( command==lengthCmd )
//	{ target->SetLength(lengthCmd->GetNewDoubleValue(newValue)); }
//}

//G4String CADPhysicsTrappingMessenger::GetCurrentValue(G4UIcommand * command)
//{
	// Print the 'trapping length' by invoking '?process/trapping/length' from the command line.
//	G4String cv;
//	if( command==lengthCmd )
//	{ cv = lengthCmd->ConvertToString(target->GetLength()); }
//	return cv;
//}
