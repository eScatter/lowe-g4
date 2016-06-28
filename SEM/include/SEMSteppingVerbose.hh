class SEMSteppingVerbose;

#ifndef SEMSteppingVerbose_h
#define SEMSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"
#include "CADPhysicsUnits.hh"

class SEMSteppingVerbose : public G4SteppingVerbose 
{
 public:
   
  SEMSteppingVerbose();
 ~SEMSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

#endif
