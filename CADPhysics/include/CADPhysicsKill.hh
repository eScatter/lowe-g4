// CADPhysicsKill.hh
//
// A very useful little process that kills a particle (currently works for electrons) after it has travelled a fixed
// distance that can be set by the user via a messenger.
// The positions where the particles die are written to a file 'killout.dat'. 

#ifndef CADPhysicsKill_h
#define CADPhysicsKill_h 1
#include "G4VContinuousDiscreteProcess.hh"
#include "G4Electron.hh"
#include <sstream>
#include "globals.hh"
#include "CADPhysicsUnits.hh"

class CADPhysicsKillMessenger;

class CADPhysicsKill : public G4VContinuousDiscreteProcess

{
public:

	CADPhysicsKill(const G4String& processName="kill");

	~CADPhysicsKill();

	G4bool IsApplicable ( const G4ParticleDefinition& );
	// returns true for electrons 

	void PrintInfoDefinition();
	// Invoked by BuildPhysicsTable().

	G4double PostStepGetPhysicalInteractionLength(
		const G4Track& track,
		G4double   previousStepSize,
		G4ForceCondition* condition
		);

	G4double GetContinuousStepLimit(const G4Track& aTrack,
		G4double previousStepSize,
		G4double currentMinimumStep,
		G4double& currentSafety); 

	G4double GetMeanFreePath(const G4Track& aTrack,
		G4double previousStepSize,
		G4ForceCondition* condition);

	G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep);

	inline void SetLength(G4double kl) {killlength=kl;}// Method (called by the messenger) to set killlength
	inline G4double GetLength() {return killlength;}// Corresponding 'get' method

private:
	//  Assignment operator hidden as private
	CADPhysicsKill & operator = (const CADPhysicsKill &right);
	CADPhysicsKill ( const CADPhysicsKill &);
	CADPhysicsKillMessenger* messenger;

	G4double killlength;

	std::ofstream output;

};

inline G4bool CADPhysicsKill::IsApplicable(const G4ParticleDefinition& particle)
{
	return ( &particle == G4Electron::Electron() );
}

inline G4double CADPhysicsKill::PostStepGetPhysicalInteractionLength(
	const G4Track&,//track
	G4double   previousStepSize,
	G4ForceCondition* condition
	)
{
	if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
		// Beginning of tracking (or just after DoIt of this process)
		theNumberOfInteractionLengthLeft = 1;//EK, NOT random but a fixed length!
		currentInteractionLength = killlength;
	} else {
		// subtract NumberOfInteractionLengthLeft 
		SubtractNumberOfInteractionLengthLeft(previousStepSize);
		if(theNumberOfInteractionLengthLeft<0.)
			theNumberOfInteractionLengthLeft=perMillion;        
	}

	// condition is set to "Not Forced"
	*condition = NotForced;

	G4double value = theNumberOfInteractionLengthLeft * currentInteractionLength;
	return value;
}

#endif
