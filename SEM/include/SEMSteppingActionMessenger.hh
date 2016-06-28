// SEMSteppingActionMessenger
//
#ifndef SEMSteppingActionMessenger_h
#define SEMSteppingActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SEMSteppingAction;
class G4UIcmdWithABool;

class SEMSteppingActionMessenger: public G4UImessenger
{
  public:
    SEMSteppingActionMessenger(SEMSteppingAction*);
    ~SEMSteppingActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMSteppingAction* target;

  private: //commands
    G4UIcmdWithABool*  checkCmd;

};

#endif


