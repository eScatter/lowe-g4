// Jul 31 2007: Added dealing with depth information (minimum z-coordinate)
// Jul 27 2007: Different way of dealing with multiple hits (EK)
// May 21 2007: Added new version of angular histogram including energy selection

#include "SEMPlane.hh"
#include "SEMPlaneMessenger.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "SEMDetectorConstruction.hh"
//#include "SEMPrimaryGeneratorAction.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "SEMTrackingAction.hh"
#include "CADPhysicsUserTrackInfo.hh"

SEMPlane::SEMPlane(G4String name) :SEMGeneralDetector(name)
{
  G4String HCname;
  HCname = name+"Collection";
  collectionName.insert(HCname);
  SEMGeneralDetector::SetDetectorType(fPlane);
  HCID = -1;
  PlaneTuple.Clear();
  CountEHPairs = false;
  messenger = new SEMPlaneMessenger(this);
  directional = false; // Default the detector is non-directional;
  direction = G4ThreeVector(0.0,0.0,1.0); // Default the detector is sensitive for particles moving in +z direction
  AllowMultipleHits = false; // Default the detector detects only the first hit for every particle
}

SEMPlane::~SEMPlane(){
  delete messenger;
}

void SEMPlane::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new SEMPlaneHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool SEMPlane::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
#ifdef SEMPlaneVerbose
  G4cerr << "Planar detector hit" << G4endl;
#endif
  G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  //if (charge !=0) return true; // TEMPORARY HACK FOR HANS MULDER TEST
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  //EB: This is not such a good idea, needs to be sorted below
  //if (CountEHPairs) preStepPoint = aStep->GetPostStepPoint();// EK, this will make particles get detected at the actual
  // position where they were killed in a detector, rather than at the starting position of the last step.
  G4ThreeVector worldPos;
  G4ThreeVector momentumDirection;
  G4TouchableHistory* theTouchable;

  // EB: CountEHPairs is turned on by the option detectStopped. In that case we want to report the end position of the
  //     particle, therefore we use the postStepPoint in that case.
  if (CountEHPairs) {
     worldPos = postStepPoint->GetPosition();
     momentumDirection = postStepPoint->GetMomentumDirection();
     theTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());
  } else {
     worldPos = preStepPoint->GetPosition();
     momentumDirection = preStepPoint->GetMomentumDirection();
     theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  }

  G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
  G4ThreeVector localMom = theTouchable->GetHistory()->GetTopTransform().TransformAxis(momentumDirection);

//  Example of how to compute the local surface normal.
//  - Find the normal in the local coordinate system of the solid
//  G4ThreeVector localNormal = (G4ThreeVector)(theTouchable->GetSolid()->SurfaceNormal(localPos));
//  G4cerr << "Normal = " << localNormal << G4endl;
//  - Transform it to the global coordinate system
//  G4ThreeVector Normal = theTouchable->GetHistory()->GetTopTransform().TransformAxis(localNormal);
//  - Compute the inner product (should always be negative as we already know we are moving into the detector)
//  G4double detinproduct = Normal.x()*momentumDirection.x()+Normal.y()*momentumDirection.y() +Normal.z()*momentumDirection.z();
//
//  G4cout << "Position  = " << worldPos.x()/nanometer << " " << worldPos.y()/nanometer << " " << worldPos.z()/nanometer << G4endl;
//  G4cout << "Normal    = " << Normal.x() << " " << Normal.y() << " " << Normal.z() << G4endl;
//  G4cout << "Direction = " << momentumDirection.x() << " " << momentumDirection.y() << " " << momentumDirection.z() << G4endl;
//  G4cout << detinproduct << G4endl;

  // It's only a hit if we are moving up into the detector 
  if (directional) {
     G4double innprod = 0;
     innprod = momentumDirection.x()*direction.x() + momentumDirection.y()*direction.y() + momentumDirection.z()*direction.z();
     if (innprod < 0.0) return true;
//        && (momentumDirection.z() < 0.)) return true;
  }

  G4double kinenergy;

  // Store the kinetic energy. If we are doing detectStopped this will give zero kinetic energy as the particle gets killed.
  if (CountEHPairs) {
     kinenergy = postStepPoint->GetKineticEnergy();
  } else {
     kinenergy = preStepPoint->GetKineticEnergy();
  }

#ifdef SEMPlaneVerbose
  G4VPhysicalVolume* theMotherPhysical = theTouchable->GetVolume(1); // mother
  G4int copyNo = theMotherPhysical->GetCopyNo();
  G4cout << "Hit detected in plane detector with copyNo = " << copyNo << G4endl; 
#endif

  // Determine the potential at this location and compute escape energy from sample
  static G4RunManager* run = G4RunManager::GetRunManager();
  SEMDetectorConstruction* detectorConstruction = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();
  G4double V;
  V = 0.0;
  // Add the inner potential
  G4Material* Material;
  if (CountEHPairs) {
     Material = postStepPoint->GetMaterial();
  } else {
     Material = preStepPoint->GetMaterial();
  }
  G4double workFunction = 0;
  G4double affinity = 0;
  G4double fermiEnergy = 0;
  if (Material) {
     // Only do something if GetMaterial returned a nonzero result. 
     // A NULL result happens if the particle stepped out of the world volume and we are using the post step point.
     G4MaterialPropertiesTable* aMPT = Material->GetMaterialPropertiesTable();

     if (aMPT) {
        if(aMPT->ConstPropertyExists("WORKFUNCTION")) workFunction = aMPT->GetConstProperty("WORKFUNCTION");
        if(aMPT->ConstPropertyExists("AFFINITY")) affinity = aMPT->GetConstProperty("AFFINITY");
        if(aMPT->ConstPropertyExists("FERMIENERGY")) fermiEnergy = aMPT->GetConstProperty("FERMIENERGY");
        if(!workFunction) workFunction = affinity;
     }
  }
  V += workFunction + fermiEnergy;

  SEMPlaneHit* aHit = new SEMPlaneHit();
  aHit->SetWorldPos(worldPos);
  aHit->SetLocalPos(localPos);
  aHit->SetMomentumDirection(momentumDirection);
  aHit->SetKinEnergy(kinenergy);
  aHit->SetTotalEnergy(kinenergy+charge*V);
  aHit->SetDepositedEnergy(aStep->GetTotalEnergyDeposit());
  if (CountEHPairs) { 
     aHit->SetTime(postStepPoint->GetGlobalTime());
  } else {
     aHit->SetTime(preStepPoint->GetGlobalTime());
  }
  aHit->SetLVAtVertex(aStep->GetTrack()->GetLogicalVolumeAtVertex());
  aHit->SetVertexPos(aStep->GetTrack()->GetVertexPosition());
  if (!(postStepPoint->GetPhysicalVolume())) { // Particle has left the world...
     aHit->SetTrackID(-1*aStep->GetTrack()->GetTrackID()); // Signal that particle is leaving the detector
  } else {
     if (preStepPoint->GetPhysicalVolume()->GetName() == postStepPoint->GetPhysicalVolume()->GetName() ) {
        aHit->SetTrackID(aStep->GetTrack()->GetTrackID());
     } else {
        aHit->SetTrackID(-1*aStep->GetTrack()->GetTrackID()); // Signal that particle is leaving the detector
     }
  }
  aHit->SetParentID(aStep->GetTrack()->GetParentID());
  CADPhysicsUserTrackInfo *trackinfo = (CADPhysicsUserTrackInfo*) aStep->GetTrack()->GetUserInformation();
  aHit->SetDepth(trackinfo->GetMinZ());
  aHit->SetRho(trackinfo->GetMaxRho());

  G4String p;

// Store the creator process, PE if it is one of the original primaries.
 if (aStep->GetTrack()->GetCreatorProcess()!=NULL) {
    aHit->SetCreatorProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
 } else {
    if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "e-" ) {
       p="PE";
    } else {
       p="Primary";
    }
    aHit->SetCreatorProcess(p);
 }

// Classify based on energy.
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "e-" ) {
    if (kinenergy<50.0*eV) {
       p="SE";
    } else {
       p="BSE";
    }
  } else {
    p=aStep->GetTrack()->GetDefinition()->GetParticleName();
  }

// Classify based on history. 
//  if (aStep->GetTrack()->GetTrackID()==1) {
//     p="BSE";
//  } else if (aStep->GetTrack()->GetParentID()==1) {
//     p="SE1"; // Label the particles that are directly generated by the primary beam as SE1 (this is not the literature definition!)
//  } else {
//     p="SE2"; // All others are SE2
//  }
  aHit->SetCategory(p);

  hitsCollection->insert(aHit);
  
  return true;
}

void SEMPlane::EndOfEvent(G4HCofThisEvent* HCE)
{
   SEMPlaneHitsCollection* Plane = 0;
   if (HCE) {
      Plane = (SEMPlaneHitsCollection*)(HCE->GetHC(HCID));
   }

   if (Plane) {
      static G4RunManager* run = G4RunManager::GetRunManager();
      const G4Event* evt = run->GetCurrentEvent();

      // Get the starting position of the primary particle resulting in this detector hit
      //SEMPrimaryGeneratorAction* primgen = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();
      //G4ThreeVector pos = primgen->GetPosition();

      G4int oldTID=-99999;
      int m_hit = Plane->entries();
      int i1;
      for (int i2=0;i2<m_hit;i2++) {
         if (CountEHPairs) {
            i1=m_hit-1-i2; // Walk through the list in reverse order in order to store the last point in the list of hits
         } else {
            i1=i2; // The first hit in the list will do...
         }
         SEMTupleElement newelement;
         SEMPlaneHit* bHit = (*Plane)[i1];
         G4bool  discard;
         G4int   TID;
         discard=false;
         //EK/EB: This routine is called once per event, hence EID is the same for all hits, we only have to compare track ID's
         //We assume that all hits the HitsCollection caused by one particle are located next to each other in the list.
         //This is true as each particle track is first completed before starting a new track.
         //Therefore we can simply discard a hit if the previous hit had the same trackid-eventid combination.
         //In fact, we're doing this after each event so we only need to check trackid.
         TID = bHit->GetTrackID();

         // We signal that the particle leaves the detector using a negative TrackID. For counting electron-hole pair
         // generation we want to discard those.
         if ((CountEHPairs) && (TID<0)) discard = true;
         TID = abs(TID);
         if (TID==oldTID) discard = true; // This prevents counting the same particle again
         oldTID=TID;

         if ((!discard) || (AllowMultipleHits)) {
            newelement.toten = bHit->GetTotalEnergy();
            newelement.kinen = bHit->GetKinEnergy();
            newelement.depositen = bHit->GetDepositedEnergy();
            newelement.momentumdirection = bHit->GetMomentumDirection();
            newelement.localposition = bHit->GetLocalPos();
            newelement.worldposition = bHit->GetWorldPos();
            newelement.vertexposition = bHit->GetVertexPos();
            newelement.lv =  bHit->GetLVAtVertex();
            newelement.process = bHit->GetCreatorProcess();
            newelement.category = bHit->GetCategory();
            //if (newelement.category=="SE1" || newelement.category=="SE2") {
			if (newelement.category=="SE") {
               fSE++;
            } else if (newelement.category=="BSE") {
               fBSE++;
            }
            newelement.ParentID = bHit->GetParentID();
            newelement.TrackID = TID;
            newelement.EventID = evt->GetEventID();
            newelement.Time = bHit->GetTime();
            newelement.Depth = bHit->GetDepth();
            newelement.Rho = bHit->GetRho();
            PlaneTuple.Add(newelement);
            //G4cerr << "Plane    hit: " << newelement.EventID << " " << newelement.TrackID << G4endl;

         }
      }
   }
}

void SEMPlane::OutputCounts(string filename){
   SEMGeneralDetector::OutputCounts(filename,PlaneTuple,fPE);
}

void SEMPlane::OutputHits(string filename){
   SEMGeneralDetector::OutputHits(filename,PlaneTuple);
}

void SEMPlane::OutputDeposits(string filename,G4double mindeposit){
   SEMGeneralDetector::OutputDeposits(filename,PlaneTuple,mindeposit);
}

void SEMPlane::OutputEnergyHistogram(string filename,G4double binsize){
   binsize *= eV;
   SEMGeneralDetector::OutputEnergyHistogram(filename,PlaneTuple,binsize,fPE);
}

void SEMPlane::OutputEnergyPDF(string filename){
   PlaneTuple.SortEnergy();
   SEMGeneralDetector::OutputEnergyPDF(filename,PlaneTuple,fPE);
}
void SEMPlane::OutputGammaPDF(string filename){
   PlaneTuple.SortEnergy();
   SEMGeneralDetector::OutputGammaPDF(filename,PlaneTuple,fPE);
}


void SEMPlane::OutputRhoHistogram(string filename,G4double binsize){
   binsize *= nanometer;
   SEMGeneralDetector::OutputRhoHistogram(filename,PlaneTuple,binsize,fPE);
}
void SEMPlane::OutputRhoPDF(string filename){
   PlaneTuple.SortRho();
   SEMGeneralDetector::OutputRhoPDF(filename,PlaneTuple,fPE);
}

void SEMPlane::OutputAngleHistogram(string filename,G4ThreeVector normal,G4int nbin){
   SEMGeneralDetector::OutputAngleHistogram(filename,PlaneTuple,fPE,normal,nbin);
}
