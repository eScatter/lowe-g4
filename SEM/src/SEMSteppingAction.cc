// 31-07-2007: EB Added depth and radius statistics (min z, max r)

#include "SEMSteppingAction.hh"
#include "SEMSteppingActionMessenger.hh"


#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"

#include "SEMDetectorConstruction.hh"
#include "SEMPrimaryGeneratorAction.hh"
#include "SEMGeneralDetector.hh"

#include "CADPhysicsUserTrackInfo.hh"


//Uncomment this line if you want to perform a check on the raytracing using the energy of the particle
//#define DoCheck

SEMSteppingAction::SEMSteppingAction()
{ 
   lastTrackID = -1; 
   lastdiff = 0.0; 
   stepnr = 0;

   messenger = new SEMSteppingActionMessenger(this);

}
SEMSteppingAction::~SEMSteppingAction()
{ 
   delete messenger;
}

void SEMSteppingAction::UserSteppingAction(const G4Step* theStep)
{ 
   // This code insures that a particle is killed if it hits BlackHole material 

   G4Track* theTrack = theStep->GetTrack();
   G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
   G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
   G4LogicalVolume* thePreLV = thePrePV->GetLogicalVolume();
   G4Material* thePreMAT = thePreLV->GetMaterial();


   if (thePreMAT->GetName()=="BlackHole"){
      theTrack->SetTrackStatus(fStopAndKill);
   }
   if (thePreMAT->GetName()=="BlackHole2") {
      // Kill particles that are almost undisturbed and will not reach the polepieces 
      G4ThreeVector p = thePrePoint->GetMomentumDirection();
      if (p.z() > 0.98) theTrack->SetTrackStatus(fStopAndKill);
   }

   // Record the deposited energy for all detectors...
   if (thePreLV->GetSensitiveDetector() != 0) {
      SEMGeneralDetector* detector = (SEMGeneralDetector *) thePreLV->GetSensitiveDetector();
      detector->AddEDeposit(theStep->GetTotalEnergyDeposit());
   }

   // Store lowest z-coordinate and largest r-coordinate (always w.r.t. z-axis)
   CADPhysicsUserTrackInfo *trackinfo = (CADPhysicsUserTrackInfo*) theTrack->GetUserInformation();
   if (theTrack->GetCurrentStepNumber()==0) {
      G4double x = theStep->GetPreStepPoint()->GetPosition().x();
      G4double y = theStep->GetPreStepPoint()->GetPosition().y();
      G4double z = theStep->GetPreStepPoint()->GetPosition().z();
      G4double rho = x*x+y*y;
      trackinfo->SetMinZ(z);
      trackinfo->SetMaxRho(rho);
   } else {
      G4double x = theStep->GetPostStepPoint()->GetPosition().x();
      G4double y = theStep->GetPostStepPoint()->GetPosition().y();
      G4double z = theStep->GetPostStepPoint()->GetPosition().z();
      G4double rho = x*x+y*y;
      if (z < trackinfo->GetMinZ() ) trackinfo->SetMinZ(z);
      if (rho > trackinfo->GetMaxRho() ) trackinfo->SetMaxRho(rho);
   }
}

void SEMSteppingAction::CompareEnergy(const G4Step* theStep) {

   G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
   G4double charge = theStep->GetTrack()->GetDefinition()->GetPDGCharge();
   G4int TrackID = theStep->GetTrack()->GetTrackID();

   G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
   if ((evt==0) && (lastTrackID < 0) ){
      lastTrackID = TrackID;
      output.open("energy_check.dat");
      output << "Step#    Position                Particle         Change   Change in      Potential    Kinetic  " << G4endl;
      output << "Step#    (mm,mm,mm)               energy  (eV)     (eV)    this step (eV) energy (eV)  energy (eV)" << G4endl;
   } else if ((evt ==0) && (TrackID == lastTrackID)) {
 
      output.open("energy_check.dat",ios::app);
   } else {
      G4cerr << "Only one particle allowed in field check." << G4endl;
      G4RunManager::GetRunManager()->AbortRun(false);
      G4cerr << "Run aborted" << G4endl;
      return;
   }

   double pos[3],gpos[3];

   // Start energy at the position of the gun
   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMDetectorConstruction*  Detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

   SEMPrimaryGeneratorAction* primgen = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();

   G4ThreeVector GunPosition = primgen->GetPosition();
   G4double      GunEnergy = primgen->GetParticleEnergy();

   gpos[0] = GunPosition.x();
   gpos[1] = GunPosition.y();
   gpos[2] = GunPosition.z();


   G4double StartEnergy = GunEnergy;

   // Present energy. The kinetic energy is obtained during raytracing by integration.
   //                 The FEM package computes the potential. In computing the E-field the gradient of
   //                 the potential is taken, thus introducing a numerical error. By checking energy conservation
   //                 we know the severity of this numerical error. The output of this routine can then be
   //                 used as a guideline for improving the FEM grid.

   G4ThreeVector Position = thePrePoint->GetPosition();
   G4double KineticEnergy = thePrePoint->GetKineticEnergy();

   pos[0] = Position.x();
   pos[1] = Position.y();
   pos[2] = Position.z();


   G4double CurrentEnergy = KineticEnergy;

   output << stepnr++ << scientific << setprecision(4) << " " 
          << setw(11) << pos[0]/mm << " " << setw(11) << pos[1]/mm << " " << setw(11) << pos[2]/mm << " " 
          << setw(11) << CurrentEnergy/eV << " " 
          << setw(11) << StartEnergy/volt - CurrentEnergy/eV << " " 
          << setw(11) << lastdiff - (StartEnergy/volt - CurrentEnergy/eV) << " "
          << setw(11) << KineticEnergy/eV << G4endl;
   
   lastdiff=StartEnergy/volt-CurrentEnergy/eV;
   output.close();
}

