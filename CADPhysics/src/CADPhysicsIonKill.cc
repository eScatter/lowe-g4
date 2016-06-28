// CADPhysicsIonKill.cc
//
#include "G4ios.hh"
#include "CADPhysicsIonKill.hh"

CADPhysicsIonKill::CADPhysicsIonKill(const G4String& processName,
									 G4ProcessType type)
									 : G4VDiscreteProcess(processName, type)
{
	if ( verboseLevel > 0) {
		G4cout << GetProcessName() << " is created " << G4endl;
	}
}

CADPhysicsIonKill::~CADPhysicsIonKill(){}

G4VParticleChange* CADPhysicsIonKill::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
	aParticleChange.Initialize(aTrack);
	aParticleChange.ProposeTrackStatus(fStopAndKill);
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


