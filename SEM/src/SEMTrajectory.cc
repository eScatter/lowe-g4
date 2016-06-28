#include "SEMTrajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

G4Allocator<SEMTrajectory> myTrajectoryAllocator;

SEMTrajectory::SEMTrajectory()
{
   fpParticleDefinition = 0;
   ParticleName = "";
   PDGCharge = 0;
   PDGEncoding = 0;
   fTrackID = 0;
   fParentID = 0;
   positionRecord = 0;
   logenergy = 1;
   startenergy = 1.*eV;
   newpolyenergy = 0.5*eV;
}

SEMTrajectory::SEMTrajectory(const G4Track* aTrack)
{
   fpParticleDefinition = aTrack->GetDefinition();
   ParticleName = fpParticleDefinition->GetParticleName();
   PDGCharge = fpParticleDefinition->GetPDGCharge();
   PDGEncoding = fpParticleDefinition->GetPDGEncoding();
   fTrackID = aTrack->GetTrackID();
   fParentID = aTrack->GetParentID();
   positionRecord = new SEMTrajectoryPointContainer();
   positionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
   startenergy = 1.*keV;
   G4double energy = aTrack->GetKineticEnergy();
   logenergy = G4int(floor(log(startenergy/energy)/log(2.)));
   newpolyenergy = startenergy*pow(0.5,logenergy);
}

SEMTrajectory::SEMTrajectory(SEMTrajectory & right)
    : G4VTrajectory()
{
  ParticleName = right.ParticleName;
  fpParticleDefinition = right.fpParticleDefinition;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  logenergy = right.logenergy;
  startenergy = right.startenergy;
  newpolyenergy = right.newpolyenergy;
  positionRecord = new SEMTrajectoryPointContainer();
  for(int i=0;i<(int)right.positionRecord->size();i++)
  {
    G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.positionRecord))[i]);
    positionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
  }
}

SEMTrajectory::~SEMTrajectory()
{
  size_t i;
  for(i=0;i<positionRecord->size();i++){
    delete  (*positionRecord)[i];
  }
  positionRecord->clear();

  delete positionRecord;
}

void SEMTrajectory::ShowTrajectory() const
{
   G4cout << G4endl << "TrackID =" << fTrackID
        << ":ParentID=" << fParentID << G4endl;
   G4cout << "Particle name : " << ParticleName
        << "  Charge : " << PDGCharge << G4endl;
   G4cout << "  Current trajectory has " << positionRecord->size()
        << " points." << G4endl;

   for( size_t i=0 ; i < positionRecord->size() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
       G4cout << "Point[" << i << "]"
            << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   }
}

void SEMTrajectory::ShowTrajectory(std::ostream& o) const
{
    G4VTrajectory::ShowTrajectory(o);
}


void SEMTrajectory::DrawTrajectory(G4int /*i_mode*/) const
{

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   G4Polyline pPolyline;
   for (int i = 0; i < (int)positionRecord->size() ; i++) {
     G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
     pos = aTrajectoryPoint->GetPosition();
     pPolyline.push_back( pos );
   }

   G4Colour colour(0.,0.,0.);    // Pitch black
   if(fpParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(0.,1.,1.);     // Cyan
   else if(fpParticleDefinition==G4Electron::ElectronDefinition())
   {
	  //colour = G4Colour(1.,0.,0.);    // Red
	  //if (smallenergy) colour = G4Colour(0.,0.,1.); //Blue?
      if (logenergy>0)
	  {
		 G4double blueamount = 1 - 1/(1+0.5*logenergy);
	     colour = G4Colour(1.-blueamount,0.,blueamount);
	  } else {
		 G4double yellowamount = -1.*logenergy/4;
		 if(1<yellowamount) yellowamount = 1;
         colour = G4Colour(1.,yellowamount,0);
	  }
   }
   // EK: other options removed
   //G4cout << "Drawing electron... \n" ;
   G4VisAttributes attribs(colour);
   pPolyline.SetVisAttributes(attribs);
   if(pVVisManager) pVVisManager->Draw(pPolyline);
}

void SEMTrajectory::AppendStep(const G4Step* aStep)
{
   G4TrajectoryPoint* point = new G4TrajectoryPoint(aStep->GetPostStepPoint()->
								GetPosition());
   G4double energy = aStep->GetPostStepPoint()->GetKineticEnergy();
   G4int newlogenergy = logenergy;
   if(0!=energy) {
	   newlogenergy = (int)(floor(log(startenergy/energy)/log(2.)));
   }
   if (logenergy!=newlogenergy) {
		  positionRecord->push_back( point );
		  DrawTrajectory();
		  positionRecord->clear();
		  logenergy=newlogenergy;
   }
   positionRecord->push_back( point );
}

G4ParticleDefinition* SEMTrajectory::GetParticleDefinition()
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void SEMTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  SEMTrajectory* seco = (SEMTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  {
    positionRecord->push_back((*(seco->positionRecord))[i]);
  }
  delete (*seco->positionRecord)[0];
  seco->positionRecord->clear();

}


