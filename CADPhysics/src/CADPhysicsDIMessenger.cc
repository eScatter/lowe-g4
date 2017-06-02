// CADPhysicsDIMessenger.cc
//

#include "CADPhysicsDIMessenger.hh"
#include "CADPhysicsDI.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

CADPhysicsDIMessenger::CADPhysicsDIMessenger(CADPhysicsDI * mpga)
:target(mpga)
{
	DIDir = new G4UIdirectory("/process/di/");
	DIDir->SetGuidance("Controls for the CADPhysicsDI process.");

	secondariesCmd = new G4UIcmdWithABool("/process/di/secondaries",this);
	secondariesCmd->SetGuidance("Generate secondaries in DI.");
	secondariesCmd->SetParameterName("flg",true);
	secondariesCmd->SetDefaultValue(true);

	xrayCmd = new G4UIcmdWithABool("/process/di/xray",this);
	xrayCmd->SetGuidance("Generate X-rays (through atomic deexcitation) in DI.");
	xrayCmd->SetParameterName("flg",true);
	xrayCmd->SetDefaultValue(true);

	augerCmd = new G4UIcmdWithABool("/process/di/auger",this);
	augerCmd->SetGuidance("Generate Auger electrons in DI.");
	augerCmd->SetParameterName("flg",true);
	augerCmd->SetDefaultValue(true);

	rangecutCmd = new G4UIcmdWithABool("/process/di/rangecut",this);
	rangecutCmd->SetGuidance("Switch for application of the range cut in DI.");
	rangecutCmd->SetParameterName("flg",true);
	rangecutCmd->SetDefaultValue(true);

	verboseCmd = new G4UIcmdWithAnInteger("/process/di/verbose",this);
	verboseCmd->SetGuidance("Verbosity level of CADPhysicsDI.");
	verboseCmd->SetGuidance("0 - silent except for essential warning messages");
	verboseCmd->SetGuidance("1 - (default) basic information during initialization");
	verboseCmd->SetGuidance("2 - more extensive information during initialization");
	verboseCmd->SetGuidance("3 - also print inverse mean free path tables to screen during initialization");
	verboseCmd->SetGuidance("4 and higher - additionally with (debugging) output during simulations");
	verboseCmd->SetParameterName("level",true);
	verboseCmd->SetRange("level>=0");
	verboseCmd->SetDefaultValue(1);

	energycutCmd = new G4UIcmdWithADoubleAndUnit("/process/di/energycut",this);
	energycutCmd->SetGuidance("Lowest allowable energy (rel. to vacuum) for an electron in a liquid/solid.");
	energycutCmd->SetGuidance("Negative values allow electrons to be tracked even if they have too low energy");
	energycutCmd->SetGuidance("to escape into vacuum.");
	energycutCmd->SetParameterName("e",true);
	energycutCmd->SetDefaultValue(0.);
	energycutCmd->SetDefaultUnit("MeV");

	resetcounterCmd = new G4UIcommand("/process/di/resetcounter",this);
	resetcounterCmd->SetGuidance("Print the number of electron/hole pairs generated so far and reset the counter.");

  outputCmd = new G4UIcmdWithABool("/process/di/output",this);
  outputCmd->SetGuidance("Generate the outputfile for the trapped electrons.");
  outputCmd->SetParameterName("flg",false);
  outputCmd->SetDefaultValue(false);

}

CADPhysicsDIMessenger::~CADPhysicsDIMessenger()
{
	delete secondariesCmd;
	delete xrayCmd;
	delete augerCmd;
	delete rangecutCmd;
	delete verboseCmd;
	delete DIDir;
	delete energycutCmd;
	delete resetcounterCmd;
  delete outputCmd;
}

void CADPhysicsDIMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
	if( command==secondariesCmd )
	{ target->SetGenerateSecondaries(secondariesCmd->GetNewBoolValue(newValue)); }
	if( command==xrayCmd )
	{ target->SetGenerateXrays(xrayCmd->GetNewBoolValue(newValue)); }
	if( command==augerCmd )
	{ target->SetGenerateAugers(augerCmd->GetNewBoolValue(newValue)); }
	if( command==rangecutCmd )
	{ target->SetRangeCut(rangecutCmd->GetNewBoolValue(newValue)); }
	if( command==verboseCmd )
	{ target->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); }
	if( command==energycutCmd )
	{ target->SetEnergyCut(energycutCmd->GetNewDoubleValue(newValue)); }
	if( command==resetcounterCmd )
	{ target->ResetCounter();}
  if( command==outputCmd )
	{ target->SetOutput(outputCmd->GetNewBoolValue(newValue)); }

}

G4String CADPhysicsDIMessenger::GetCurrentValue(G4UIcommand * command)
{
	// Uses the 'get' methods of CADPhysicsDI to print the current values of the given parameters.
	// These can be accessed by replacing the first '/' by a question mark; e.g. '?process/di/verbose' to
	// request the current verbosity level of the process.
	G4String cv;
	if( command==secondariesCmd )
	{ cv = secondariesCmd->ConvertToString(target->GetGenerateSecondaries()); }
	if( command==xrayCmd )
	{ cv = xrayCmd->ConvertToString(target->GetGenerateXrays()); }
	if( command==augerCmd )
	{ cv = augerCmd->ConvertToString(target->GetGenerateAugers()); }
	if( command==rangecutCmd )
	{ cv = rangecutCmd->ConvertToString(target->GetRangeCut()); }
	if( command==verboseCmd )
	{ cv = verboseCmd->ConvertToString(target->GetVerboseLevel()); }
	if( command==energycutCmd )
	{ cv = energycutCmd->ConvertToString(target->GetEnergyCut()); }

	return cv;
}
