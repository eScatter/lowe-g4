// 31-07-2007: EB added depth information 

#include "SEMTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "SEMTrajectory.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "SEMRunAction.hh"

#include "CADPhysicsUserTrackInfo.hh"

SEMTrackingAction::SEMTrackingAction()
{
}

SEMTrackingAction::~SEMTrackingAction()
{;}

void SEMTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
   //EK, commented out:
   //fpTrackingManager->SetStoreTrajectory(true);

   //EB: Original code
   G4bool storetrajectory = fpTrackingManager->GetStoreTrajectory();
   if(storetrajectory) fpTrackingManager->SetTrajectory(new SEMTrajectory(aTrack));

   //EB: Added code for finding maximum depth of track
   CADPhysicsUserTrackInfo* trackinfo = 0;
   if (aTrack->GetUserInformation()) {
	   trackinfo = (CADPhysicsUserTrackInfo*) aTrack->GetUserInformation();
   } else {
	   trackinfo = new CADPhysicsUserTrackInfo;
	   fpTrackingManager->SetUserTrackInformation(trackinfo);
   }
   trackinfo->SetMinZ(aTrack->GetPosition().z());
   G4double x,y;
   x=aTrack->GetPosition().x();
   y=aTrack->GetPosition().y();
   trackinfo->SetMaxRho(x*x+y*y);
}


void SEMTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
   if (aTrack->GetTrackID() == 1) {
      // Code for mapping the ends of all primary tracks. 
      // If /scan/mapOn is active, this code returns the colour of the lv where the track ended.
      SEMRunAction* ra = (SEMRunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
      if (ra->GetMapTrackEnd()) {
         if (! (aTrack->GetVolume()->GetLogicalVolume()->GetVisAttributes()) ) {
            ra->SetMapColour(G4Colour(0.0,0.0,0.0,0.0)); // Volume has no visual attributes, make it black
         } else {
            ra->SetMapColour(aTrack->GetVolume()->GetLogicalVolume()->GetVisAttributes()->GetColour());
         }
      }
   }
}
