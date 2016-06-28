// SEMGeneralDetectorMessenger
//
#ifndef SEMGeneralDetectorMessenger_h
#define SEMGeneralDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SEMGeneralDetector;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class SEMGeneralDetectorMessenger: public G4UImessenger
{
  public:
    SEMGeneralDetectorMessenger(SEMGeneralDetector* mpga);
    ~SEMGeneralDetectorMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMGeneralDetector* target;

  private: //commands
    G4UIcmdWithAnInteger*  ModelCmd;
    G4UIcmdWithABool*  EwindowCmd;
    G4UIcmdWithADoubleAndUnit*  minECmd;
    G4UIcmdWithADoubleAndUnit*  maxECmd;
    G4UIcmdWithABool*  AwindowCmd;
    G4UIcmdWithADoubleAndUnit*  minACmd;
    G4UIcmdWithADoubleAndUnit*  maxACmd;
    G4UIcmdWith3Vector*         dirAVecCmd;
    G4UIcmdWithAString*         selectOutputCmd;
    G4UIcmdWithAString*         deselectOutputCmd;

};

#endif


