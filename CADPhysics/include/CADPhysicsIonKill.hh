// CADPhysicsIonKill.hh
//
// Kills an ion if it is not in a gas.

#ifndef CADPhysicsIonKill_h
#define CADPhysicsIonKill_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

#include "G4Material.hh" 
#include "G4MaterialCutsCouple.hh"

class CADPhysicsIonKill : public G4VDiscreteProcess 
{

public:
	CADPhysicsIonKill(const G4String& processName = "IonKill",
		G4ProcessType type = fOptical);
	~CADPhysicsIonKill();

public:

	G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
	// Returns true -> 'is applicable' only for an ion

	G4double PostStepGetPhysicalInteractionLength(
		const G4Track& track,
		G4double   previousStepSize,
		G4ForceCondition* condition
		);

	G4double GetMeanFreePath(const G4Track& ,
		G4double ,
		G4ForceCondition* condition);
	// Returns infinity; i. e. the process does not limit the step,
	// but sets the 'Forced' condition for the DoIt to be invoked at
	// every step.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
		const G4Step&  aStep);
	// This is the method implementing the boundary processes.

private:

};

inline
G4bool CADPhysicsIonKill::IsApplicable(const G4ParticleDefinition& 
									   aParticleType)
{
	return ( aParticleType.GetParticleType() == "nucleus");
}

inline G4double CADPhysicsIonKill::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,G4double   , G4ForceCondition* condition)
{
	*condition = NotForced;
	G4double value = DBL_MAX;
	const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
	const G4Material* material = couple->GetMaterial();
	if(material->GetState()!=kStateGas) {
		// Set the physical interaction length to zero if the ion is not in a gas - this has the consequence that the ion will be killed
		// at the end of this step.
		value=0.;
	}
	return value;
}

inline G4double CADPhysicsIonKill::GetMeanFreePath(const G4Track& ,G4double ,G4ForceCondition* condition)
// This method is not even called by PostStepGetPhysicalInteractionLength but it returns a force condition just to be sure.
{
	*condition = Forced;
	return DBL_MAX;
}

#endif /* CADPhysicsIonKill_h */
