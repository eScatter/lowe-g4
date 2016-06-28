// CADPhysicsUserTrackInfo.cc
//

#include "CADPhysicsUserTrackInfo.hh"

CADPhysicsUserTrackInfo::CADPhysicsUserTrackInfo(G4double ee)
: excitationenergy(ee),// The excitation energy can be set but just calling the constructor. Convenient for CADPhysicsGasScattering.
minZ(DBL_MAX),
maxRho(-1.0),
uservar1(0.),
uservar2(0.),
uservar3(0.)
{}

CADPhysicsUserTrackInfo::~CADPhysicsUserTrackInfo() {}
