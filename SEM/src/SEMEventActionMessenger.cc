// SEMEventActionMessenger.cc
//

#include "SEMEventActionMessenger.hh"

#include "SEMEventAction.hh"
#include "G4UIcmdWithAnInteger.hh"

SEMEventActionMessenger::SEMEventActionMessenger(SEMEventAction * mpga)
:target(mpga)
{
  verboseCmd = new G4UIcmdWithAnInteger("/detectors/verbose",this);
  verboseCmd->SetGuidance("Verbose level for each event.");
  verboseCmd->SetGuidance(" Event summary will be displayed for every 'level' events.");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetRange("level>=0");
  verboseCmd->SetDefaultValue(1);

  randomCmd = new G4UIcmdWithAnInteger("/random/seed",this);
  randomCmd->SetGuidance("Set random seed.");
  randomCmd->SetGuidance("Default is to use the system clock for the initial seed.");
  randomCmd->SetParameterName("seed",true);
  randomCmd->SetRange("seed>=0");
}

SEMEventActionMessenger::~SEMEventActionMessenger()
{
  delete verboseCmd;
  delete randomCmd;
}

void SEMEventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==verboseCmd ) { target->SetVerbose(verboseCmd->GetNewIntValue(newValue)); }
  if( command==randomCmd ) { target->SetSeed(randomCmd->GetNewIntValue(newValue)); }
}

G4String SEMEventActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==verboseCmd ) { cv = verboseCmd->ConvertToString(target->GetVerbose()); }
  if( command==randomCmd ) { cv = randomCmd->ConvertToString(target->GetSeed()); }

  return cv;
}
