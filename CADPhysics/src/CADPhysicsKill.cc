// CADPhysicsKill.cc
//

#include "CADPhysicsKill.hh"
#include "CADPhysicsKillMessenger.hh"

using namespace std;

CADPhysicsKill::CADPhysicsKill(const G4String& processName)
     : G4VContinuousDiscreteProcess(processName)
  { 
    SetVerboseLevel(1);
	killlength = DBL_MAX;
	messenger = new CADPhysicsKillMessenger(this);
  }

CADPhysicsKill::~CADPhysicsKill()
{
	output.close();
}

G4double CADPhysicsKill::GetContinuousStepLimit(
                                   const G4Track&,//track
                                   G4double,
                                   G4double,//currentMinimumStep,
                                   G4double&)
{
  // Not implemented
  return DBL_MAX;
}

G4VParticleChange* CADPhysicsKill::AlongStepDoIt(
                                       const G4Track& track,const G4Step&)//track,step
{
  //Not implemented
  aParticleChange.Initialize(track);
  return &aParticleChange;
}

G4VParticleChange* CADPhysicsKill::PostStepDoIt(
                                               const G4Track& track,
                                               const G4Step& step)
{ 
  // First thing to do: kill the particle.
  aParticleChange.Initialize(track);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  G4ThreeVector KillPosition = step.GetPostStepPoint()->GetPosition();
  G4double x = KillPosition.x()/nanometer;
  G4double y = KillPosition.y()/nanometer;
  G4double z = KillPosition.z()/nanometer;
  G4double r = sqrt(x*x+y*y);
  if(!output.is_open()) {
    output.open ("killout.dat");
    output << setprecision(10);
	output << "length\tx\ty\tr\tz" << G4endl;
  }
  output << killlength/nanometer << "\t" << x << "\t" << y << "\t" << r << "\t" << z << G4endl;
  
  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


G4double CADPhysicsKill::GetMeanFreePath(
                                            const G4Track&,// track,
                                            G4double,
                                            G4ForceCondition*)
{
	return DBL_MAX;
}


void CADPhysicsKill::PrintInfoDefinition()
{
	G4cout << "CADPhysicsKill: Kill a particle after a fixed track length.\n";
}
