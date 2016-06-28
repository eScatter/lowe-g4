// CADPhysicsIonBoundary.cc
//

#include "CADPhysicsIonBoundary.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "CADPhysicsIonBoundaryMessenger.hh"
#include "G4Material.hh"
#include "G4TransportationManager.hh"

CADPhysicsIonBoundary::CADPhysicsIonBoundary(const G4String& processName,
											 G4ProcessType type)
											 : G4VDiscreteProcess(processName, type),
											 gamma(0.1)// Default version of gamma
{
	if ( verboseLevel > 0) {
		G4cout << GetProcessName() << " is created " << G4endl;
	}
	trackid = -1;
	settokill = false;
	messenger = new CADPhysicsIonBoundaryMessenger(this);
}

CADPhysicsIonBoundary::~CADPhysicsIonBoundary(){}

G4VParticleChange* CADPhysicsIonBoundary::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
	aParticleChange.Initialize(aTrack);

	G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
	if (settokill && pPreStepPoint->GetStepStatus() == fGeomBoundary) {// Even if the flag 'settokill' is true,
		// an extra check is done to see if the particle indeed started its current step at a boundary.
		aParticleChange.ProposeTrackStatus(fStopAndKill);// Kill the particle
		settokill = false;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
	settokill = false;

	// Next, check if the particle is currently at a boundary.
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
	if (pPostStepPoint->GetStepStatus() != fGeomBoundary)
	{
		// If not, do not do anything else
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	// Now, if the material of the volume that the particle is about to enter is a gas, we do not need to do 
	// anything either.
	Material2 = pPostStepPoint -> GetMaterial();
	if (Material2->GetState()==kStateGas)// Check the physical state of the material
	{
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	// Now we know that the particle is about to enter a non-gaseous material. 
	// Set a flag to kill the particle after the next step and store its track ID for later reference.
	settokill = true;
	trackid = aTrack.GetTrackID();

	// The following part of the code only serves to find the surface normal, and to test whether
	// it is pointing in the proper direction (which, for our purposes, means pointing out of 
	// the second volume back into the first volume). This test is done by calculating the inner product
	// of the surface normal with the current direction of the particle; this inner product should be negative.
	// Note: this piece of code was copied from G4OpBoundaryProcess.cc (Geant4 v. 8.0p01). It may be worthwhile
	// to check in future versions of Geant4 to see if this part of code has been improved (simplified).
	const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
	G4ThreeVector OldMomentum = aParticle->GetMomentumDirection();
	G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
	G4Navigator* theNavigator =	G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	G4ThreeVector theLocalPoint = theNavigator->GetGlobalToLocalTransform().TransformPoint(theGlobalPoint);
	G4ThreeVector theLocalNormal;
	G4bool valid;
	theLocalNormal = theNavigator->GetLocalExitNormal(&valid);
	if (valid) {
		theLocalNormal = -theLocalNormal;
	} else {
		G4cerr << " CADPhysicsIonBoundary/PostStepDoIt(): "
			<< " The Navigator reports that it returned an invalid normal" 
			<< G4endl;
	}
	theGlobalNormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
	if ( verboseLevel > 3 ) {
		G4cout << "theLocalPoint=" << theLocalPoint << G4endl;
		G4cout << "theGlobalPoint=" << theGlobalPoint << G4endl;
		G4cout << "theLocalNormal=" << theLocalNormal << G4endl;
		G4cout << "theGlobalNormal=" << theGlobalNormal << G4endl;
	}
	if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4DEBUG_OPTICAL
		G4cerr << " CADPhysicsIonBoundary/PostStepDoIt(): "
			<< " theGlobalNormal points the wrong direction "
			<< G4endl;
#endif
		theGlobalNormal = -theGlobalNormal;
	}
	// End of this part

	// Next, take care of creating a secondary electron.
	G4double OldEnergy = aParticle->GetKineticEnergy();
	if (G4UniformRand()<gamma && OldEnergy>150.*eV) {// Currently the secondary electron emission coefficient
		// is equal to a fixed value gamma for ion energies above 150 eV and zero below. Of course much 
		// refinement is still possible here.
		// A random number determines whether a secondary electron is actually generated.

		// Determine the direction of the secondary electron
		G4double R = G4UniformRand();
		G4double cosChidaughter = sqrt(R);// In this manner cos(Chi) has a probability distribution that is proportional to cos(Chi).
		// This is the so-called cosine distribution.
		G4double sinChidaughter = sqrt(1.-R);
		G4double Phi = twopi*G4UniformRand();// Phi has a uniform random distribution
		G4double dirxd = sinChidaughter*cos(Phi), diryd = sinChidaughter*sin(Phi), dirzd = cosChidaughter;
		G4DynamicParticle* EmittedElectron = new G4DynamicParticle();// Create the new particle

		G4double energydaughter = 1.;// The SE energy is fixed at 1 eV, for now
		EmittedElectron->SetKineticEnergy(energydaughter * eV);
		G4ThreeVector daughterdir = (G4ThreeVector(dirxd,diryd,dirzd)).unit();
		// Align the SE direction towards the surface normal
		daughterdir.rotateUz(theGlobalNormal);
		EmittedElectron->SetMomentumDirection(daughterdir);
		EmittedElectron->SetDefinition(G4Electron::Electron());
		aParticleChange.SetNumberOfSecondaries(1);
		aParticleChange.AddSecondary(EmittedElectron);
	}
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4double CADPhysicsIonBoundary::GetMeanFreePath(const G4Track& ,
												G4double ,
												G4ForceCondition* condition)
{
	*condition = Forced;
	return DBL_MAX;
}

