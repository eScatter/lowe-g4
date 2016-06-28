#ifndef SEMScanMessenger_h
#define SEMScanMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SEMRunAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class SEMScanMessenger: public G4UImessenger
{
  public:
    SEMScanMessenger(SEMRunAction* scanAct);
    ~SEMScanMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMRunAction * scanManager;

  private: //commands
    G4UIdirectory *              scanDirectory;
    G4UIcommand *                mapCmd;
    G4UIcommand *                scanCmd;
    G4UIcommand *                imageCmd;
    G4UIcommand *                scanOutputCmd;
    G4UIcommand *                imageOutputCmd;
    G4UIcommand *                imageASCIIOutputCmd;
    G4UIcmdWith3Vector *         scanXdirCmd;
    G4UIcmdWith3Vector *         scanYdirCmd;
    G4UIcmdWith3Vector *         mapRotCmd;
    G4UIcmdWith3VectorAndUnit *  scanOriginCmd;
    G4UIcmdWithAnInteger *       scanNxCmd;
    G4UIcmdWithAnInteger *       scanNyCmd;
    G4UIcmdWithADoubleAndUnit *  scanLxCmd;
    G4UIcmdWithADoubleAndUnit *  scanLyCmd;
    G4UIcmdWithAString *         prefixCmd;
};

#endif


