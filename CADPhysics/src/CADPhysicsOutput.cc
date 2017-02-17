// CADPhysicsOutput.cc
//

#include "CADPhysicsOutput.hh"

namespace CADPhysicsOutput
{

G4String CADPhysicsOutput(const G4Track& track,
                 const G4Step& step)
                 {
                   // First thing to do: stop the particle.
                   G4double kineticEnergy = track.GetKineticEnergy();
                   G4int EID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
                   G4int TID = track.GetTrackID();
                   G4double totalEnergy = track.GetTotalEnergy();
                   G4double xOrigin = track.GetVertexPosition().x();
                   G4double yOrigin = track.GetVertexPosition().y();
                   G4double zOrigin = track.GetVertexPosition().z();
                   G4double xDirection = track.GetMomentumDirection().x();
                   G4double yDirection = track.GetMomentumDirection().y();
                   G4double zDirection = track.GetMomentumDirection().z();
                   CADPhysicsUserTrackInfo *trackinfo = (CADPhysicsUserTrackInfo*) track.GetUserInformation();
                   G4double MaxDepth = trackinfo->GetMinZ();
                   G4double MaxR = sqrt(trackinfo->GetMaxRho());

                   G4ThreeVector TrapPosition = step.GetPostStepPoint()->GetPosition();
                   G4double x = TrapPosition.x();
                   G4double y = TrapPosition.y();
                   G4double z = TrapPosition.z();

                   std::stringstream output;
                   output << EID << "\t" << TID << "\t" << totalEnergy/eV << "\t" << kineticEnergy/eV << "\t"
                          << x/nanometer << "\t" << y/nanometer << "\t" << z/nanometer << "\t"
                          << xOrigin/nanometer << "\t" << yOrigin/nanometer << "\t" << zOrigin/nanometer << "\t"
                          << xDirection/nanometer << "\t" << yDirection/nanometer << "\t" << zDirection/nanometer << "\t"
                          << MaxDepth/nanometer << "\t" << MaxR/nanometer;
                   return output.str();
                 }
}
