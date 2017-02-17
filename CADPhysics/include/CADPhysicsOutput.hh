// CADPhysicsOutput.hh
//
// A class to output information of trapped or absorbed electrons

#ifndef CADPhysicsOutput_h
#define CADPhysicsOutput_h 1
#include "globals.hh"
#include <sstream>
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "CADPhysicsUserTrackInfo.hh"

// namespace CADPhysicsOutput;

namespace CADPhysicsOutput

{
	G4String CADPhysicsOutput(const G4Track& track,
                   const G4Step& step);
}

#endif
