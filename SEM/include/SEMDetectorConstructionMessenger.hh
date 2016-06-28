#ifndef SEMDetectorConstructionMessenger_h
#define SEMDetectorConstructionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SEMDetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;

class SEMDetectorConstructionMessenger: public G4UImessenger
{
  public:
    SEMDetectorConstructionMessenger(SEMDetectorConstruction* mpga);
    ~SEMDetectorConstructionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    SEMDetectorConstruction *target;

  private: //commands
    G4UIdirectory           *detectorsDirectory;
    G4UIcommand             *writeCmd;
	G4UIcommand*  materialCmd;
	G4UIcmdWithAString* gasfractionsCmd;
	G4UIcommand*  gasdensityCmd;
	G4UIcmdWithAString* getdensityCmd;
	G4UIcommand*  printtreeCmd;
};

#endif


