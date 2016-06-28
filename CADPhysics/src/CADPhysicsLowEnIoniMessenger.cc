// CADPhysicsLowEnIoniMessenger.cc
//

#include "CADPhysicsLowEnIoniMessenger.hh"
#include "CADPhysicsLowEnergyIonisation.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

CADPhysicsLowEnIoniMessenger::CADPhysicsLowEnIoniMessenger(CADPhysicsLowEnergyIonisation * mpga)
:target(mpga)
{
	LEIDir = new G4UIdirectory("/process/loweioni/");
	LEIDir->SetGuidance("Controls for the CADPhysicsLowEnergyIonisation process.");

	usecutCmd = new G4UIcmdWithABool("/process/loweioni/usecut",this);
	usecutCmd->SetGuidance("Apply the energy cut for electron generation in CADPhysicsLowEnergyIonisation.");
	usecutCmd->SetParameterName("flg",true);
	usecutCmd->SetDefaultValue(true);

	generateionsCmd = new G4UIcmdWithABool("/process/loweioni/generateions",this);
	generateionsCmd->SetGuidance("Create ions for tracking in CADPhysicsLowEnergyIonisation.");
	generateionsCmd->SetParameterName("flg",true);
	generateionsCmd->SetDefaultValue(true);

	ecutCmd = new G4UIcmdWithADoubleAndUnit("/process/loweioni/ecut",this);
	ecutCmd->SetGuidance("Set the cut energy value for electrons in CADPhysicsLowEnergyIonisation.");
	ecutCmd->SetParameterName("Energy",true,true);
	ecutCmd->SetDefaultUnit("GeV");
}

CADPhysicsLowEnIoniMessenger::~CADPhysicsLowEnIoniMessenger()
{
	delete usecutCmd;
	delete generateionsCmd;
	delete ecutCmd;
	delete LEIDir;
}

void CADPhysicsLowEnIoniMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
	if( command==usecutCmd )
	{ target->SetUseCut(usecutCmd->GetNewBoolValue(newValue)); }
	if( command==generateionsCmd )
	{ target->SetGenerateIons(generateionsCmd->GetNewBoolValue(newValue)); }
	if( command==ecutCmd )
	{ target->SetCutForLowEnSecElectrons(ecutCmd->GetNewDoubleValue(newValue)); }
}

G4String CADPhysicsLowEnIoniMessenger::GetCurrentValue(G4UIcommand * command)
{
	// Methods to print the current values of the various parameters to screen.
	// Can be executed from the command line by entering
	// '?process/loweioni/usecut',
	// '?process/loweioni/generateions', or
	// '?process/loweioni/ecut'.
	G4String cv;
	if( command==usecutCmd )
	{ cv = usecutCmd->ConvertToString(target->GetUseCut()); }
	if( command==generateionsCmd )
	{ cv = generateionsCmd->ConvertToString(target->GetGenerateIons()); }
	if( command==ecutCmd )
	{ cv = ecutCmd->ConvertToString(target->GetCutForLowEnSecElectrons()); }
	return cv;
}
