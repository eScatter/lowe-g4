#ifndef SEMSteppingAction_h
#define SEMSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <fstream>
using namespace std;

class G4Step;
class G4StepPoint;
class SEMSteppingActionMessenger;

class SEMSteppingAction : public G4UserSteppingAction
{
  public:
    SEMSteppingAction();
   ~SEMSteppingAction();

    void UserSteppingAction(const G4Step*);

    SEMSteppingActionMessenger *messenger;

  private: // Needed for field energy checking option
    ofstream output;
    G4int    stepnr;
    G4double lastdiff;
    G4int    lastTrackID;
    void CompareEnergy(const G4Step*);

};

#endif

