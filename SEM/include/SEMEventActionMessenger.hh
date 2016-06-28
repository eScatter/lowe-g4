// SEMEventActionMessenger
//
#ifndef SEMEventActionMessenger_h
#define SEMEventActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SEMEventAction;
class G4UIcmdWithAnInteger;

class SEMEventActionMessenger: public G4UImessenger
{
  public:
    SEMEventActionMessenger(SEMEventAction* mpga);
    ~SEMEventActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMEventAction* target;

  private: //commands
    G4UIcmdWithAnInteger*  verboseCmd;
    G4UIcmdWithAnInteger*  randomCmd;

};

#endif


