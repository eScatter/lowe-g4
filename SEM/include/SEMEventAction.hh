// SEMEventAction.hh
//
#ifndef SEMEventAction_h
#define SEMEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "Randomize.hh"

class SEMEventActionMessenger;


class SEMEventAction : public G4UserEventAction
{
  public:
    SEMEventAction();
    virtual ~SEMEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    SEMEventActionMessenger* messenger;
    G4int verboseLevel;

    G4int RandomSeed;
 
  public:
    inline void SetVerbose(G4int val) { verboseLevel = val; }
    inline G4int GetVerbose() const { return verboseLevel; }

    inline void SetSeed(G4int val) { RandomSeed = val; CLHEP::HepRandom::setTheSeed(RandomSeed);}
    inline G4int GetSeed() const { return RandomSeed; }
};

#endif
