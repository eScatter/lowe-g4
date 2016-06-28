#include "SEMPrimaryGeneratorAction.hh"
#include "SEMPrimaryGeneratorActionMessenger.hh"
#include "SEMDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4RunManager.hh"

SEMPrimaryGeneratorAction::SEMPrimaryGeneratorAction()
{
  particleGun = new G4GeneralParticleSource();
  particleGun2 = new G4ParticleGun(); // Used for generating particles using a file generated with OutputHits

  gpsgun=true; // Default is to use the gps gun

  // Add a /gps/ion command that also defines a transportation process on the fly
  // The "right" way to do this would be through the physicslist for every ion type separately
  messenger = new SEMPrimaryGeneratorActionMessenger(this); 
}

SEMPrimaryGeneratorAction::~SEMPrimaryGeneratorAction()
{
  delete particleGun;
  delete particleGun2;
}

void SEMPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (gpsgun) { 
     particleGun->GeneratePrimaryVertex(anEvent);
  } else {
     G4int EID = anEvent->GetEventID();
     if (EID >= (G4int) (KinEnergy.size())) {
        // We ran out of hits in the file, switch to gps source and warn the user
//        G4cout << "WARNING: No more hits in the file, switching to gps source!" << G4endl; 
//        gpsgun=true;
//        particleGun->GeneratePrimaryVertex(anEvent);
        G4cout << "WARNING: No more hits in the file, aborting run" << G4endl; 
        G4RunManager::GetRunManager()->AbortRun();

     } else {
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        particleGun2->SetParticleDefinition(particleTable->FindParticle(ParticleType[EID]));

//        if (ParticleType[EID]=="electron") {
//           particleGun->SetParticleDefinition(G4Electron::ElectronDefinition());
//           G4cerr << "Electron" << G4endl;
//        } else if (ParticleType[EID]=="gamma") {
//           particleGun->SetParticleDefinition(G4Gamma::GammaDefinition());
//           G4cerr << "Photon" << G4endl;
//        }

        particleGun2->SetParticleEnergy(KinEnergy[EID]);
        particleGun2->SetParticlePosition(Position[EID]);
        particleGun2->SetParticleMomentumDirection(Direction[EID]);
        particleGun2->GeneratePrimaryVertex(anEvent);
     }
  }
}

// The following methods are only needed for the gps gun

void SEMPrimaryGeneratorAction::SetPosition(G4ThreeVector pos)
{ 
//  particleGun->SetParticlePosition(position);  // Old version
  position = pos;
  particleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);
}

G4ThreeVector SEMPrimaryGeneratorAction::GetPosition()
{
  position = particleGun->GetParticlePosition();
  return position;
}

G4double SEMPrimaryGeneratorAction::GetParticleEnergy()
{
  return (particleGun->GetParticleEnergy());
}

G4int SEMPrimaryGeneratorAction::GetNumberOfParticles()
{
   return particleGun->GetNumberOfParticles();
}
void SEMPrimaryGeneratorAction::SetNumberOfParticles(G4int n)
{
   particleGun->SetNumberOfParticles(n);
}

G4ThreeVector SEMPrimaryGeneratorAction::GetParticleMomentumDirection()
{
   return particleGun->GetParticleMomentumDirection();
}

G4SingleParticleSource* SEMPrimaryGeneratorAction::GetCurrentSource()
{
   return particleGun->GetCurrentSource();
}

// This one is for the hits gun

void SEMPrimaryGeneratorAction::ReadHitsFile(G4String filename)
{
   G4String line;

   ifstream hitsFile;
   hitsFile.open(filename);
   if (!hitsFile) {
      G4Exception("Could not open hits file","FileOpenError",FatalException,"SEMPrimaryGeneratorAction(1)");
   }
   getline(hitsFile,line); // Get the header
   getline(hitsFile,line); // Get the header
   while (!hitsFile.eof()) {
      G4int    EID,TID,PID;
      G4double Etot, Ekin,Time,depth,rho;
      G4ThreeVector position,vposition,direction;
      G4String lv,process,category;
      getline(hitsFile,line);
      if (line !="") {
         istringstream is(line);
         is >> EID >> TID >> PID >> Etot >> Ekin >> Time >> position >> vposition >> direction >> depth >> rho >> lv >> process >> category;
         //G4cerr << setw(7) << EID << " "
         //       << setw(7) << TID << " "
         //       << setw(7) << PID << " "
         //       << scientific << setprecision(6)
         //       << setw(15) << Etot/eV
         //       << setw(15) << Ekin/eV
         //       << setw(15) << Time/ns
         //       << setw(15) << position.x()
         //       << setw(15) << position.y()
         //       << setw(15) << position.z()
         //       << setw(15) << direction.x()
         //       << setw(15) << direction.y()
         //       << setw(15) << direction.z()
         //       << G4endl;
         KinEnergy.push_back(Ekin*eV);
         Position.push_back(position*nanometer);
         Direction.push_back(direction);
         if (category=="<gamma>") {
            ParticleType.push_back("gamma");
         }else {
            ParticleType.push_back("e-");
         }

      }
   }
}


