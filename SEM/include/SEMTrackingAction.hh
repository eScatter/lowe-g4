#ifndef SEMTrackingAction_h
#define SEMTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class SEMTrackingAction : public G4UserTrackingAction
{
  public:
    SEMTrackingAction();
    virtual ~SEMTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
};

#endif
