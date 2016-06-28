// SEMSteppingActionMessenger.cc
//

#include "SEMSteppingActionMessenger.hh"

#include "SEMSteppingAction.hh"
#include "G4UIcmdWithABool.hh"

#include "G4RunManager.hh"
#include "SEMDetectorConstruction.hh"

SEMSteppingActionMessenger::SEMSteppingActionMessenger(SEMSteppingAction * mpga)
:target(mpga)
{
    static G4RunManager* run = G4RunManager::GetRunManager();
    SEMDetectorConstruction*  Detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();
	checkCmd = 0;
}

SEMSteppingActionMessenger::~SEMSteppingActionMessenger()
{
  if (checkCmd) delete checkCmd;
}

void SEMSteppingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
}

G4String SEMSteppingActionMessenger::GetCurrentValue(G4UIcommand * /*command*/)
{
  G4String cv;

  return cv;
}
