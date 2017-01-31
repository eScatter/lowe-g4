// CADPhysicsTrapping.cc
//

#include "CADPhysicsTrapping.hh"
#include "CADPhysicsTrappingMessenger.hh"

using namespace std;

CADPhysicsTrapping::CADPhysicsTrapping(const G4String& processName)
     : G4VContinuousDiscreteProcess(processName)
  {
    SetVerboseLevel(1);
	traplength = DBL_MAX;
  trap_C = 1.; // nm^-1 value for alumina according to Ganachaud
  trap_g = 0.25; // eV^-1 value for alumina according to Ganachaud
  messenger = new CADPhysicsTrappingMessenger(this);
  }

CADPhysicsTrapping::~CADPhysicsTrapping()
{
	output.close();
}

G4double CADPhysicsTrapping::GetContinuousStepLimit(
                                   const G4Track&,//track
                                   G4double,
                                   G4double,//currentMinimumStep,
                                   G4double&)
{
  // Not implemented
  return DBL_MAX;
}

//G4VParticleChange* CADPhysicsTrapping::AlongStepDoIt(
//                                       const G4Track& track,const G4Step&)//track,step
//{
  //Not implemented
//  aParticleChange.Initialize(track);
//  return &aParticleChange;
//}

G4VParticleChange* CADPhysicsTrapping::PostStepDoIt(
                                               const G4Track& track,
                                               const G4Step& step)
{
  // First thing to do: stop the particle.
  G4double kineticEnergy = track.GetKineticEnergy();
  aParticleChange.Initialize(track);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  aParticleChange.ProposeLocalEnergyDeposit(kineticEnergy);

  G4ThreeVector TrapPosition = step.GetPostStepPoint()->GetPosition();
  G4double x = TrapPosition.x()/nanometer;
  G4double y = TrapPosition.y()/nanometer;
  G4double z = TrapPosition.z()/nanometer;
  G4double r = sqrt(x*x+y*y);
  if(!output.is_open()) {
    output.open ("trapout.dat");
    output << setprecision(10);
	output << "length\tx\ty\tr\tz" << G4endl;
  }
  output << traplength/nanometer << "\t" << x << "\t" << y << "\t" << r << "\t" << z << G4endl;

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


G4double CADPhysicsTrapping::GetMeanFreePath(
                                            const G4Track& track,
                                            G4double,
                                            G4ForceCondition*)
{
  G4double kineticEnergy = track.GetKineticEnergy();
  G4double imfp_trap = trap_C*exp(-trap_g*kineticEnergy/eV); // energy dependent cross section for trapping in nm^-1
  traplength = (millimeter/nanometer)*1./imfp_trap; // mfp in mm

	return traplength;
}


void CADPhysicsTrapping::PrintInfoDefinition()
{
	G4cout << "CADPhysicsTrapping: Trap a particle with energy dependent cross section.\n";
}
