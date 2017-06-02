// CADPhysicsTrapping.hh
// AT: based on CADPhysicsKill
//
// A process to trap electrons using the model as described by Ganachaud and Mokrani (1995) and Dapor, Ciappa and Fichtner (2010)

#ifndef CADPhysicsTrapping_h
#define CADPhysicsTrapping_h 1
#include "G4VContinuousDiscreteProcess.hh"
#include "G4Electron.hh"
#include <sstream>
#include "globals.hh"
#include "CADPhysicsUnits.hh"
#include "CADPhysicsTrappingMessenger.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"

// class CADPhysicsKillMessenger;

class CADPhysicsTrapping : public G4VContinuousDiscreteProcess

{
public:

	CADPhysicsTrapping(const G4String& processName="trapping");

	~CADPhysicsTrapping();

	G4bool IsApplicable ( const G4ParticleDefinition& );
	// returns true for electrons

	void PrintInfoDefinition();
	// Invoked by BuildPhysicsTable().

//	G4double PostStepGetPhysicalInteractionLength(
//		const G4Track& track,
//		G4double   previousStepSize,
//		G4ForceCondition* condition
//		);

	G4double GetContinuousStepLimit(const G4Track& aTrack,
		G4double previousStepSize,
		G4double currentMinimumStep,
		G4double& currentSafety);

	G4double GetMeanFreePath(const G4Track& aTrack,
		G4double previousStepSize,
		G4ForceCondition* condition);

//	G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep);

	inline void SetTrap_C(G4double kl) {trap_C=1./kl;}// Method (called by the messenger) to set trap_C
	inline void SetTrap_g(G4double kl) {trap_g=1./kl;}// Method (called by the messenger) to set trap_g
  inline void SetOutput(G4bool kl) {trap_output=kl;}// Method (called by the messenger) to set trap_output
//	inline G4double GetLength() {return traplength;}// Corresponding 'get' method

private:
	//  Assignment operator hidden as private
	CADPhysicsTrapping & operator = (const CADPhysicsTrapping &right);
	CADPhysicsTrapping ( const CADPhysicsTrapping &);
	CADPhysicsTrappingMessenger* messenger;

	G4double traplength;
	G4double trap_C;
	G4double trap_g;

  G4bool trap_output;

	std::ofstream output;

};

inline G4bool CADPhysicsTrapping::IsApplicable(const G4ParticleDefinition& particle)
{
	return ( &particle == G4Electron::Electron() );
}

#endif
