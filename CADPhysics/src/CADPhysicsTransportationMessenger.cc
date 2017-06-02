// CADPhysicsTransportationMessenger.cc
//

#include "CADPhysicsTransportationMessenger.hh"
#include "CADPhysicsTransportation.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"

CADPhysicsTransportationMessenger::CADPhysicsTransportationMessenger(CADPhysicsTransportation * mpga)
:target(mpga)
{
	transportationDir = new G4UIdirectory("/process/transportation/");
	transportationDir->SetGuidance("Controls for the CADPhysicsTransportation process.");

  outputCmd = new G4UIcmdWithABool("/process/transportation/output",this);
	outputCmd->SetGuidance("Generate the outputfile for the trapped electrons.");
	outputCmd->SetParameterName("flg",false);
	outputCmd->SetDefaultValue(false);
}

CADPhysicsTransportationMessenger::~CADPhysicsTransportationMessenger()
{
	delete transportationDir;
  delete outputCmd;
}

void CADPhysicsTransportationMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==outputCmd )
	{ target->SetOutput(outputCmd->GetNewBoolValue(newValue)); }
}
