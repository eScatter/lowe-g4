#ifndef SEMPlaneMessenger_h
#define SEMPlaneMessenger_h 1

class SEMPlane;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;
class G4UIcmdWithAString;
class G4UIcommand;
//class G4UIcmdWithAnInteger;
//class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"

class SEMPlaneMessenger: public G4UImessenger
{
  public:
    SEMPlaneMessenger(SEMPlane * mpga);
    ~SEMPlaneMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMPlane* target;

  private: //commands
    G4UIdirectory*         planeDirectory;
    G4UIcmdWithABool*      multiCmd;
    G4UIcmdWithABool*      pairCmd;
    G4UIcmdWithABool*      dirCmd;
    G4UIcmdWith3Vector*    dirVecCmd;
    G4UIcmdWithAString*    planeCountsOutputCmd;
    G4UIcmdWithAString*    planeHitsOutputCmd;
    G4UIcommand*           planeHistOutputCmd;
    G4UIcmdWithAString*    planeEnergyPDFOutputCmd;
    G4UIcmdWithAString*    planeGammaPDFOutputCmd;
    G4UIcommand*           planeRhoHistOutputCmd;
    G4UIcmdWithAString*    planeRhoPDFOutputCmd;
    G4UIcommand*           planeAngleHistOutputCmd;
    G4UIcommand*           planeDepositsOutputCmd;
    G4UIcmdWithADoubleAndUnit* dirThetaCmd;
    G4UIcmdWithADoubleAndUnit* dirPhiCmd;
    G4double dirTheta,dirPhi;

};

#endif


