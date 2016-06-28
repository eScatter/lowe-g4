// CADPhysicsIonBoundary.hh, loosely based on G4OpBoundaryProcess.hh
//
// This process takes care of killing ions when they hit surfaces of non-gaseous materials.
// Also, it can create a secondary electron that is emitted from the surface into the gas. In other words,
// it does NOT simulate the detailed ion-material interaction, but generates secondary electrons
// following a predefined distribution.
//
// Since this process only acts if the particle is at a surface, its structure is similar to that of CADPhysicseBoundary.
//

#ifndef CADPhysicsIonBoundary_h
#define CADPhysicsIonBoundary_h 1
#include "globals.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Electron.hh"
#include "CADPhysicsUnits.hh"

class G4Material;
class CADPhysicsIonBoundaryMessenger;

class CADPhysicsIonBoundary : public G4VDiscreteProcess 
{

public:
	CADPhysicsIonBoundary(const G4String& processName = "IonBoundary",
		G4ProcessType type = fOptical);
	~CADPhysicsIonBoundary();

public:
	G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
	// Returns true for particle types 'nucleus' (covers 'ordinary' ions) and 'molecule' (for the CADPhysicsWaterIon).

	G4double PostStepGetPhysicalInteractionLength(
		const G4Track& track,
		G4double   previousStepSize,
		G4ForceCondition* condition
		);
	// Special implementation that does not actually limit the step length, but ensures that PostStepDoIt is called for every step. 

	G4double GetMeanFreePath(const G4Track& ,
		G4double ,
		G4ForceCondition* condition);
	// Does the same as PostStepGetPhysicalInteractionLength. In fact, since this method is NOT used by PostStepGetPhysicalInteractionLength
	// (as in 'normal' processes), it is actually redundant. 

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
		const G4Step&  aStep);
	// This method first checks if the particle is at a surface (of a non-gaseous material), then if it is, handles creation of secondary electrons.
	// It also sets a flag to kill the ion. The ion is actually killed only in the *next* step, so that it is allowed to set one step inside the 
	// second volume. If the second volume is defined as a detector, this allows that detector to generate a hit from the ion.

	inline void SetGamma(G4double g) {gamma = g;}
	inline G4double GetGamma() {return gamma;}

private:

	CADPhysicsIonBoundaryMessenger* messenger;

	G4ThreeVector theGlobalNormal;
	G4ThreeVector theFacetNormal;

	G4Material* Material2;// The material that the ion is entering

	G4double gamma;// The ion-induced secondary electron emission coefficient
	G4bool settokill;// Flag to kill the ion in the next step
	G4int trackid;// ID of the ion track. This number is used to check if in the next call of PostStepDoIt, we are still dealing with the same
	// ion. (Normally, this should always be the case.)

};

inline G4bool CADPhysicsIonBoundary::IsApplicable(const G4ParticleDefinition& aParticleType)
{
	return ( aParticleType.GetParticleType() == "nucleus" || aParticleType.GetParticleType() == "molecule" );
}

inline G4double CADPhysicsIonBoundary::PostStepGetPhysicalInteractionLength(
	const G4Track& track,
	G4double,// previousStepSize
	G4ForceCondition* condition
	)
{
	*condition = Forced;// Force PostStepDoIt to be executed after every step
	if (settokill && trackid == track.GetTrackID() ) {
		return DBL_MIN;// Return a tiny step size if the ion was set to be killed after this step.
	}
	settokill = false;
	return DBL_MAX;// Otherwise, do not limit the step size.
}

#endif
