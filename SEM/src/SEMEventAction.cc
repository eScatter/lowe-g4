// SEMEventAction.cc
//
#include "SEMEventAction.hh"
#include "SEMEventActionMessenger.hh"

#include "G4Event.hh"
#include <time.h>

SEMEventAction::SEMEventAction()
{
  verboseLevel = 1;
  messenger = new SEMEventActionMessenger(this);

  // Initialize random number generator using system clock
  SetSeed((G4int) time(NULL));
  G4cout << "Random number generator seeded with " << GetSeed() << G4endl;
}

SEMEventAction::~SEMEventAction()
{
   delete messenger;
}

void SEMEventAction::BeginOfEventAction(const G4Event*)
{
}

void SEMEventAction::EndOfEventAction(const G4Event* evt)
{
   if (verboseLevel!=0 && evt->GetEventID() % verboseLevel == 0) {
      G4cout  << "Event " << evt->GetEventID() << G4endl;
   }
} 
