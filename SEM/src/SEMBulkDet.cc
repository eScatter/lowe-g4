#include "SEMBulkDet.hh"
#include "SEMBulkDetMessenger.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "SEMDetectorConstruction.hh"
//#include "SEMPrimaryGeneratorAction.hh"

#include "CADPhysicsUserTrackInfo.hh"

SEMBulkDet::SEMBulkDet(G4String name) :SEMGeneralDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="BulkDetCollection");
  SEMGeneralDetector::SetDetectorType(fBulk);
  HCID = -1;
  DetTuple.Clear();
  messenger = new SEMBulkDetMessenger(this);
  detectAtEnd = false; // Default the detector records the data at the startpoint of each hit
  multipleHits = false; // Default the detector detects only the first or last hit for every particle (depending on detectAtEnd)
  detectTotal = false; // Default the detector does not only give the totals for the event
  energyCut = 0.0; // Default is to store all hits that result in a positive energy deposit
}

SEMBulkDet::~SEMBulkDet(){
   delete messenger;
}

void SEMBulkDet::Initialize(G4HCofThisEvent*HCE)
{
   hitsCollection = new SEMBulkDetHitsCollection(SensitiveDetectorName,collectionName[0]);
   if ( HCID<0 ) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); 
   HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool SEMBulkDet::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
  // Flags:
  // - detectAtEnd: true if the user wants to know where particles end up. Information is given at the end point of the step

  // Works for both electrons and gamma
  G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();

  // If detectAtEnd we use the information at the end of the step
  G4StepPoint* theStepPoint;
  if (detectAtEnd) {
     theStepPoint = aStep->GetPostStepPoint();
  } else {
     theStepPoint = aStep->GetPreStepPoint();
  }

  G4ThreeVector       worldPos = theStepPoint->GetPosition();
  G4ThreeVector       momentumDirection = theStepPoint->GetMomentumDirection();
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(theStepPoint->GetTouchable());
  G4double            kinenergy =  theStepPoint->GetKineticEnergy(); 
  G4Material*         Material = theStepPoint->GetMaterial();

  // Local position
  G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  // Determine the potential at this location and compute escape energy from sample
  static G4RunManager* run = G4RunManager::GetRunManager();
  SEMDetectorConstruction* detectorConstruction = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();
  G4double V;
  V = 0.0;
  // Add the inner potential
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

  SEMBulkDetHit* aHit = new SEMBulkDetHit();

  aHit->SetWorldPos(worldPos);
  aHit->SetLocalPos(localPos);
  aHit->SetMomentumDirection(momentumDirection);
  aHit->SetKinEnergy(kinenergy);  
  aHit->SetTotalEnergy(kinenergy+charge*V);  
  aHit->SetDepositedEnergy(aStep->GetTotalEnergyDeposit());

  aHit->SetTime(theStepPoint->GetGlobalTime());

  CADPhysicsUserTrackInfo *trackinfo = (CADPhysicsUserTrackInfo*) aStep->GetTrack()->GetUserInformation();
  aHit->SetDepth(trackinfo->GetMinZ());
  aHit->SetRho(trackinfo->GetMaxRho());

  aHit->SetLVAtVertex(aStep->GetTrack()->GetLogicalVolumeAtVertex());
  aHit->SetVertexPos(aStep->GetTrack()->GetVertexPosition());
  aHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  aHit->SetParentID(aStep->GetTrack()->GetParentID());

// Store the creator process, PE if it is one of the original primaries.
 G4String p;
 if (aStep->GetPostStepPoint()->GetProcessDefinedStep()!=NULL) {
    aHit->SetCreatorProcess(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
 } else {
    p="UserLimit";
    aHit->SetCreatorProcess(p);
 }
// Store the creator process, PE if it is one of the original primaries.
// G4String p;
// if (aStep->GetTrack()->GetCreatorProcess()!=NULL) {
//    aHit->SetCreatorProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
// } else {
//    if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "e-" ) {
//       p="PE";
//    } else {
//       p="Primary";
//    }
//    aHit->SetCreatorProcess(p);

// Classify based on kinetic energy.
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "e-" ) {
    if (kinenergy<50.0*eV) {
       p="SE";
    } else {
       p="BSE";
    }
  } else {
    p=aStep->GetTrack()->GetDefinition()->GetParticleName();
  }
  aHit->SetCategory(p); 

  hitsCollection->insert(aHit);

  return true;
}

void SEMBulkDet::EndOfEvent(G4HCofThisEvent* HCE)
{
   SEMBulkDetHitsCollection* DHC1 = 0;
   if (HCE) {
      DHC1 = (SEMBulkDetHitsCollection*)(HCE->GetHC(HCID));
   }

   if (DHC1) {
      static G4RunManager* run = G4RunManager::GetRunManager();
      const G4Event* evt = run->GetCurrentEvent();

      // Get the starting position and direction of the primary particle resulting in this detector hit
      //SEMPrimaryGeneratorAction* primgen = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();
      //G4ThreeVector pos = primgen->GetPosition();
      //G4ThreeVector dir = primgen->GetParticleMomentumDirection();

      // The detector was hit, store the required statistics of this hit 
      // Flags:
      // - multipleHits  true:  every single hit is stored, the deposited energy during the step is given
      //                 false: only the first or last hit is stored, depending on the value of detectAtEnd
      //                        the total deposited energy along this track is given
      //                       
      int n_hit = DHC1->entries();
      if (n_hit !=0 ) {
         SEMTupleElement newelement; 
         SEMBulkDetHit* aHit;
         if (detectTotal) {
            // We only need totals: store the first hit and add the deposits
            aHit = (*DHC1)[0]; // Should be the first hit for this event
            newelement.toten = aHit->GetTotalEnergy();
            newelement.kinen = aHit->GetKinEnergy();
            newelement.depositen = 0.0;
            newelement.momentumdirection = aHit->GetMomentumDirection();
            newelement.localposition = aHit->GetLocalPos();
            newelement.worldposition = aHit->GetWorldPos();
            newelement.vertexposition = aHit->GetVertexPos();
            newelement.lv =  aHit->GetLVAtVertex();
            newelement.process = aHit->GetCreatorProcess();
            newelement.category = aHit->GetCategory();
            newelement.ParentID = aHit->GetParentID();
            newelement.TrackID = aHit->GetTrackID();
            newelement.EventID = evt->GetEventID();
            newelement.Time = aHit->GetTime();
            newelement.Depth = aHit->GetDepth();
            newelement.Rho = aHit->GetRho();
            for (int i=0;i<n_hit;i++) {
               aHit = (*DHC1)[i];
               newelement.depositen += aHit->GetDepositedEnergy();
            }
            DetTuple.Add(newelement);
         } else {
            G4int    oldTID=-1;
            G4bool   sameParticle=false;
            G4int    i1;
            G4double SmallEnergyDepositsFromPreviousSteps=0.0; // Stores small energy losses (< energyCut) which will be skipped. The energy is added to the next step (only for MultipleHits=true)
            // We have to walk through the list and store depending on multipleHits and detectAtEnd
            for (int i2=0;i2<n_hit;i2++) {
               if ((detectAtEnd)&&(!multipleHits)) {
                  // Walk through the list in descending order
                  i1=n_hit-1-i2; 
               } else {
                  // Walk through the list in ascending order
                  i1=i2; 
               }
               aHit = (*DHC1)[i1];

               G4int TID = aHit->GetTrackID();
               // 
               sameParticle = (TID==oldTID);
               oldTID=TID;
   
               if (multipleHits) {
                  // Every hit that resulted in a deposited energy larger than energyCut is stored
                  newelement.toten = aHit->GetTotalEnergy();
                  newelement.kinen = aHit->GetKinEnergy();
                  newelement.depositen = aHit->GetDepositedEnergy()+SmallEnergyDepositsFromPreviousSteps;
                  newelement.momentumdirection = aHit->GetMomentumDirection();
                  newelement.localposition = aHit->GetLocalPos();
                  newelement.worldposition = aHit->GetWorldPos();
                  newelement.vertexposition = aHit->GetVertexPos();
                  newelement.lv =  aHit->GetLVAtVertex();
                  newelement.process = aHit->GetCreatorProcess();
                  newelement.category = aHit->GetCategory();
                  newelement.ParentID = aHit->GetParentID();
                  newelement.TrackID = aHit->GetTrackID();
                  newelement.EventID = evt->GetEventID();
                  newelement.Time = aHit->GetTime();
                  newelement.Depth = aHit->GetDepth();
                  newelement.Rho = aHit->GetRho();
                  if (newelement.depositen < energyCut) {
                     SmallEnergyDepositsFromPreviousSteps = newelement.depositen; // Store this small deposit so we don't have to bother about this step
                  } else {
                     SmallEnergyDepositsFromPreviousSteps=0.0; 
                     DetTuple.Add(newelement);
                  }
               } else if (!sameParticle) {
                  // Only first or last hit is stored, we found a new particle so time to store the hit
                  if (i2 !=0) DetTuple.Add(newelement); // Only store the hit if we have anything to store
                  // Store the first hit for the new particle
                  newelement.toten = aHit->GetTotalEnergy();
                  newelement.kinen = aHit->GetKinEnergy();
                  newelement.depositen = aHit->GetDepositedEnergy();
                  newelement.momentumdirection = aHit->GetMomentumDirection();
                  newelement.localposition = aHit->GetLocalPos();
                  newelement.worldposition = aHit->GetWorldPos();
                  newelement.vertexposition = aHit->GetVertexPos();
                  newelement.lv =  aHit->GetLVAtVertex();
                  newelement.process = aHit->GetCreatorProcess();
                  newelement.category = aHit->GetCategory();
                  newelement.ParentID = aHit->GetParentID();
                  newelement.TrackID = aHit->GetTrackID();
                  newelement.EventID = evt->GetEventID();
                  newelement.Time = aHit->GetTime();
                  newelement.Depth = aHit->GetDepth();
                  newelement.Rho = aHit->GetRho();
               } else {
                  newelement.depositen += aHit->GetDepositedEnergy();
               }
            }  // for loop
            if (!multipleHits) DetTuple.Add(newelement); // Store the last particle too, note that the if (n_hit!=0) prevents this happening in case of no hits
         } // if (detectTotal)
      }
   }
}

void SEMBulkDet::OutputCounts(string filename){
   SEMGeneralDetector::OutputCounts(filename,DetTuple,fPE);
}
void SEMBulkDet::OutputHits(string filename){
   SEMGeneralDetector::OutputHits(filename,DetTuple);
}

void SEMBulkDet::OutputDeposits(string filename){
   G4double mindeposit = 0.0; // For this detector we must count all deposits
   SEMGeneralDetector::OutputDeposits(filename,DetTuple,mindeposit);
}

void SEMBulkDet::OutputEnergyHistogram(string filename,G4double binsize){
   binsize *= eV;
   SEMGeneralDetector::OutputEnergyHistogram(filename,DetTuple,binsize,fPE);
}
void SEMBulkDet::OutputTimeHistogram(string filename,G4double binsize){
   binsize *= ns;
   SEMGeneralDetector::OutputTimeHistogram(filename,DetTuple,binsize,fPE);
}
void SEMBulkDet::OutputEnergyPDF(string filename){
   DetTuple.SortEnergy();
   SEMGeneralDetector::OutputEnergyPDF(filename,DetTuple,fPE);
}
void SEMBulkDet::OutputRhoPDF(string filename){
   DetTuple.SortRho();
   SEMGeneralDetector::OutputRhoPDF(filename,DetTuple,fPE);
}
void SEMBulkDet::OutputAngleHistogram(string filename,G4ThreeVector normal,G4int nbin){
   SEMGeneralDetector::OutputAngleHistogram(filename,DetTuple,fPE,normal,nbin);
}

