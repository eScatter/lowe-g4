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

	trap_C_Cmd = new G4UIcmdWithADoubleAndUnit("/process/trapping/invC",this);
	trap_C_Cmd->SetGuidance("Set 1/trap_C used to calculate the trapping mean free path.");
	trap_C_Cmd->SetParameterName("invC",true);
  trap_C_Cmd->SetDefaultValue(1.); // C = 1. nm^-1 for alumina according to Ganachaud
	trap_C_Cmd->SetDefaultUnit("nanometer");
	trap_C_Cmd->SetUnitCandidates("nanometer micron mm cm m km");
	trap_C_Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  trap_g_Cmd = new G4UIcmdWithADoubleAndUnit("/process/trapping/invG",this);
	trap_g_Cmd->SetGuidance("Set 1/trap_g used to calculate the trapping mean free path.");
	trap_g_Cmd->SetParameterName("invG",true);
	trap_g_Cmd->SetDefaultValue(4.); // g = 0.25 eV^-1 for alumina according to Ganachaud
	trap_g_Cmd->SetDefaultUnit("eV");
	trap_g_Cmd->SetUnitCandidates("eV keV MeV GeV TeV");
	trap_g_Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

CADPhysicsTrappingMessenger::~CADPhysicsTrappingMessenger()
{
	delete trap_C_Cmd;
  delete trap_g_Cmd;
	delete trappingDir;
}

void CADPhysicsTrappingMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
	if( command==trap_C_Cmd )
	{ target->SetTrap_C(trap_C_Cmd->GetNewDoubleValue(newValue)); }
  if( command==trap_g_Cmd )
	{ target->SetTrap_g(trap_g_Cmd->GetNewDoubleValue(newValue)); }
}

//G4String CADPhysicsTrappingMessenger::GetCurrentValue(G4UIcommand * command)
//{
	// Print the 'trapping length' by invoking '?process/trapping/length' from the command line.
//	G4String cv;
//	if( command==lengthCmd )
//	{ cv = lengthCmd->ConvertToString(target->GetLength()); }
//	return cv;
//}
