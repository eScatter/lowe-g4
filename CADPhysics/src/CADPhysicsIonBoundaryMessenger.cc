// CADPhysicsIonBoundaryMessenger.cc
//

#include "CADPhysicsIonBoundaryMessenger.hh"
#include "CADPhysicsIonBoundary.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"

CADPhysicsIonBoundaryMessenger::CADPhysicsIonBoundaryMessenger(CADPhysicsIonBoundary * mpga)
:target(mpga)
{
	IBDir = new G4UIdirectory("/process/ionboundary/");
	IBDir->SetGuidance("Controls for the CADPhysicsIonBoundary process.");

	gammaCmd = new G4UIcmdWithADouble("/process/ionboundary/gamma",this);
	gammaCmd->SetGuidance("Set the ion-induced secondary emission coefficient.");
	gammaCmd->SetParameterName("gamma",true);
	gammaCmd->SetDefaultValue(0.);
}

CADPhysicsIonBoundaryMessenger::~CADPhysicsIonBoundaryMessenger()
{
	delete IBDir;
	delete gammaCmd;
}

void CADPhysicsIonBoundaryMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
	if( command==gammaCmd )
	{ target->SetGamma(gammaCmd->GetNewDoubleValue(newValue)); }
}

G4String CADPhysicsIonBoundaryMessenger::GetCurrentValue(G4UIcommand * command)
{
	// Print the current value of gamma (call by entering '?process/ionboundary/gamma' from the commmand line)
	G4String cv;
	if( command==gammaCmd )
	{ cv = gammaCmd->ConvertToString(target->GetGamma()); }
	return cv;
}
