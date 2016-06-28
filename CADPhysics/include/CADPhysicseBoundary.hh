// CADPhysicseBoundary.hh, loosely based on G4OpBoundaryProcess.hh
//
// Takes care of the behaviour of an electron at a volume boundary.
// The electron may get either transmitted or absorbed/reflected.
//

#ifndef CADPhysicseBoundary_h
#define CADPhysicseBoundary_h 1
#include "globals.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Electron.hh"
#include "CADPhysicsUnits.hh"

class G4Material;

enum CADPhysicseBoundaryStatus {
	Undefined,
	Refraction, 
	TotalInternalReflection,
	NotAtBoundary,
	SameMaterial, StepTooSmall,
	ParticleKilled
};

class CADPhysicseBoundary : public G4VDiscreteProcess 
{

public:
	CADPhysicseBoundary(const G4String& processName = "eBoundary",
		G4ProcessType type = fOptical);
	~CADPhysicseBoundary();

	G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
	// Returns true -> 'is applicable' only for an electron.

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
	// This method checks if the particle is at a boundary and executes the boundary processes if it is.

	CADPhysicseBoundaryStatus GetStatus() const;
	// Returns the current status.

private:

	void DielectricDielectric();
	// Private method that performs the actual check for reflection/transmission

	// Data members

	G4ThreeVector OldMomentum;
	G4double	  OldEnergy;

	G4ThreeVector NewMomentum;
	G4double	  NewEnergy;

	G4ThreeVector theGlobalNormal;
	G4ThreeVector theFacetNormal;

	G4Material* Material1;// The materials on both sides of the surface
	G4Material* Material2;

	G4double workFunction1;// Parameters of both materials
	G4double workFunction2;
	G4double affinity1;
	G4double affinity2;
	G4double fermiEnergy1;
	G4double fermiEnergy2;
	G4double bandbending1;
	G4double bandbending2;
	G4double deltaphi1;
	G4double deltaphi2;
	G4double U1;
	G4double U2;

	G4double deltaU;
	G4double Ugain;
	G4bool	 vacuuminterface;
	G4int    lasttrackID;

	CADPhysicseBoundaryStatus theStatus;

};

inline G4bool CADPhysicseBoundary::IsApplicable(const G4ParticleDefinition& 
												aParticleType)
{
	return ( &aParticleType == G4Electron::Electron() );
}

inline G4double CADPhysicseBoundary::PostStepGetPhysicalInteractionLength(
	const G4Track&,// track
	G4double,// previousStepSize
	G4ForceCondition* condition
	)
{
	*condition = Forced;
	return DBL_MAX;
}

inline CADPhysicseBoundaryStatus CADPhysicseBoundary::GetStatus() const
{
	return theStatus;
}

#endif
