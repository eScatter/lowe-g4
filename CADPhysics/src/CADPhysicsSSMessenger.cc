// CADPhysicsSSMessenger.cc
//

#include "CADPhysicsSSMessenger.hh"
#include "CADPhysicsSingleScattering.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

CADPhysicsSSMessenger::CADPhysicsSSMessenger(CADPhysicsSingleScattering * mpga)
:target(mpga)
{
	SSDir = new G4UIdirectory("/process/ss/");
	SSDir->SetGuidance("Controls for the CADPhysicsSingleScattering process.");

	multistepCmd = new G4UIcmdWithABool("/process/ss/multistep",this);
	multistepCmd->SetGuidance("Do 'multistep' software acceleration in CADPhysicsSingleScattering.");
	multistepCmd->SetParameterName("flg",true);
	multistepCmd->SetDefaultValue(true);

	diffusionstepCmd = new G4UIcmdWithABool("/process/ss/diffusionstep",this);
	diffusionstepCmd->SetGuidance("Do 'diffusionstep' software acceleration in CADPhysicsSingleScattering.");
	diffusionstepCmd->SetParameterName("flg",true);
	diffusionstepCmd->SetDefaultValue(true);

	verboseCmd = new G4UIcmdWithAnInteger("/process/ss/verbose",this);
	verboseCmd->SetGuidance("Verbosity level of CADPhysicsSingleScattering.");
	verboseCmd->SetGuidance("0 - nothing");
	verboseCmd->SetGuidance("1 - (default) some basic information during initialization");
	verboseCmd->SetGuidance("2 - more detailed information during initialization");
	verboseCmd->SetGuidance("3 or 4 - additionally, print elastic mean free paths during initialization");
	verboseCmd->SetGuidance("5 and higher - also print the mean free path for each step during simulations");
	verboseCmd->SetParameterName("level",true);
	verboseCmd->SetRange("level>=0");
	verboseCmd->SetDefaultValue(1);
}

CADPhysicsSSMessenger::~CADPhysicsSSMessenger()
{
	delete multistepCmd;
	delete diffusionstepCmd;
	delete verboseCmd;
	delete SSDir;
}

void CADPhysicsSSMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
	if( command==multistepCmd )
	{ target->SetDoMultistep(multistepCmd->GetNewBoolValue(newValue)); }
	if( command==diffusionstepCmd )
	{ target->SetDoDiffusionstep(diffusionstepCmd->GetNewBoolValue(newValue)); }
	if( command==verboseCmd )
	{ target->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); }
}

G4String CADPhysicsSSMessenger::GetCurrentValue(G4UIcommand * command)
{
	// Uses the 'get' methods of CADPhysicsSingleScattering to print the current values of the given parameters.
	// These can be accessed by replacing the first '/' by a question mark; e.g. '?process/ss/verbose' to 
	// request the current verbosity level of the process.
	G4String cv;
	if( command==multistepCmd )
	{ cv = multistepCmd->ConvertToString(target->GetDoMultistep()); }
	if( command==multistepCmd )
	{ cv = diffusionstepCmd->ConvertToString(target->GetDoDiffusionstep()); }
	if( command==verboseCmd )
	{ cv = verboseCmd->ConvertToString(target->GetVerboseLevel()); }
	return cv;
}
