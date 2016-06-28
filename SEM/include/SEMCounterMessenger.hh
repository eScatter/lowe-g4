#ifndef SEMCounterMessenger_h
#define SEMCounterMessenger_h 1

class SEMCounter;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcommand;

#include "G4UImessenger.hh"

class SEMCounterMessenger: public G4UImessenger
{
  public:
    SEMCounterMessenger(SEMCounter * mpga);
    ~SEMCounterMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMCounter* target;

  private: //commands
    G4UIdirectory*         counterDirectory;
    G4UIcmdWithAString*    counterCountsOutputCmd;
    G4UIcmdWithAString*    counterHitsOutputCmd;
    G4UIcommand*           counterHistOutputCmd;
    G4UIcmdWithAString*    counterEnergyPDFOutputCmd;
    G4UIcmdWithAnInteger*  counterVerboseCmd;

};

#endif


