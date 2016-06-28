#ifndef SEMBulkDetMessenger_h
#define SEMBulkDetMessenger_h 1

class SEMBulkDet;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcommand;

#include "G4UImessenger.hh"

class SEMBulkDetMessenger: public G4UImessenger
{
  public:
    SEMBulkDetMessenger(SEMBulkDet * mpga);
    ~SEMBulkDetMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMBulkDet* target;

  private: //commands
    G4UIdirectory*         detDirectory;
    G4UIcmdWithABool*      detMultipleHitsCmd;
    G4UIcmdWithABool*      detDetectAtEndCmd;
    G4UIcmdWithABool*      detDetectTotalCmd;
    G4UIcmdWithADoubleAndUnit* detEnergyCutCmd;
    G4UIcmdWithAString*    detCountsOutputCmd;
    G4UIcmdWithAString*    detHitsOutputCmd;
    G4UIcommand*           detHistOutputCmd;
    G4UIcmdWithAString*    detEnergyPDFOutputCmd;
    G4UIcmdWithAString*    detRhoPDFOutputCmd;
    G4UIcommand*           detAngleHistOutputCmd;
    G4UIcommand*           detTimeHistOutputCmd;
    G4UIcommand*           detDepositsOutputCmd;

};

#endif


