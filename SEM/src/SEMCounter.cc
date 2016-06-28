// Oct 16 2008: Started Counter type detector based on Plane type detector
//              The idea is that this detector only records the deposited energy in every event

// Jul 31 2007: Added dealing with depth information (minimum z-coordinate)
// Jul 27 2007: Different way of dealing with multiple hits (EK)
// May 21 2007: Added new version of angular histogram including energy selection

#include "SEMCounter.hh"
#include "SEMCounterMessenger.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "SEMDetectorConstruction.hh"
#include "SEMPrimaryGeneratorAction.hh"

#include "SEMTrackingAction.hh"
#include "CADPhysicsUserTrackInfo.hh"

SEMCounter::SEMCounter(G4String name) :SEMGeneralDetector(name)
{
  G4String HCname;
  HCname = name+"Collection";
  collectionName.insert(HCname);
  SEMGeneralDetector::SetDetectorType(fCounter);
  HCID = -1;
  CounterTuple.Clear();
  nhit=0;
  EHit=0.0;
  DepositedTillNow=0.0;
  verboseLevel=0;
  messenger = new SEMCounterMessenger(this);
}

SEMCounter::~SEMCounter(){
  delete messenger;
}

void SEMCounter::Initialize(G4HCofThisEvent* /*HCE*/)
{
  // EB: this line was used for checking the energy balance within the code. It is removed as we should count the cumulative deposit for a RUN.
  //  SetEDeposit(0.0); // Clear deposited energy before each event;
  EHit=0.0; // Clear energy of hit.
  DepositedTillNow=GetEDeposit();
  lastTID=-1;
}

G4bool SEMCounter::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
  if (GetLVname() != aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()) {
     // This particle came from outside
     G4int thisTID = aStep->GetTrack()->GetTrackID();
     if (lastTID != thisTID) {
        EHit += aStep->GetPreStepPoint()->GetKineticEnergy();
        lastTID=thisTID;
     }
  } else if (GetLVname() != aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() ) {
    // This particle originated in the detector but is leaving subtract its kinetic energy
     G4int thisTID = aStep->GetTrack()->GetTrackID();
     if (lastTID != thisTID) {
        EHit -= aStep->GetPostStepPoint()->GetKineticEnergy();
        lastTID=thisTID;
     }
  }
  return true;
}

void SEMCounter::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
   if (lastTID > 0) {
      static G4RunManager* run = G4RunManager::GetRunManager();
      const G4Event* evt = run->GetCurrentEvent();
      SEMTupleElement newelement;

      newelement.toten = GetEDeposit()-DepositedTillNow; // Total deposited energy 
      newelement.kinen = EHit; // Energy of incoming particle(s)
      newelement.momentumdirection = G4ThreeVector(0.0,0.0,0.0);
      newelement.localposition = G4ThreeVector(0.0,0.0,0.0);
      newelement.worldposition = G4ThreeVector(0.0,0.0,0.0);
      newelement.vertexposition = G4ThreeVector(0.0,0.0,0.0);
      newelement.lv =  "NA";
      newelement.process = "NA";
      newelement.category = "NA";
      newelement.ParentID = 0;
      newelement.TrackID = 0;
      newelement.EventID = evt->GetEventID();
      newelement.Time = 0.0;
      newelement.Depth = 0.0;
      newelement.Rho = 0.0;

      CounterTuple.Add(newelement);
      nhit++;
      if (verboseLevel>0) {
         G4cout << "Hit# " << nhit << " Kinetic energy of incoming particle(s) " << EHit/eV << " eV" << G4endl;
         G4cout << "Energy deposited = " << (GetEDeposit()-DepositedTillNow)/eV << " eV " << G4endl;
      }
   }
}

void SEMCounter::OutputCounts(string filename){
   SEMGeneralDetector::OutputCounts(filename,CounterTuple,fPE);
}

void SEMCounter::OutputHits(string filename){
   SEMGeneralDetector::OutputHits(filename,CounterTuple);
}

void SEMCounter::OutputEnergyHistogram(string filename,G4double binsize){
   binsize *= eV;
   SEMGeneralDetector::OutputEnergyHistogram(filename,CounterTuple,binsize,fPE);
}

void SEMCounter::OutputEnergyPDF(string filename){
   CounterTuple.SortEnergy();
   SEMGeneralDetector::OutputEnergyPDF(filename,CounterTuple,fPE);
}
