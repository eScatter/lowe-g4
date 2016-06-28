class SEMTrajectory;

#ifndef SEMTrajectory_h
#define SEMTrajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"                 // Include from 'system'
#include <vector>             //
#include <iostream>             //
#include "globals.hh"               // Include from 'global'
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "G4TrajectoryPoint.hh"     // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Colour.hh"
#include "CADPhysicsUnits.hh"

typedef std::vector<G4VTrajectoryPoint*> SEMTrajectoryPointContainer;

class G4Polyline;                   // Forward declaration.

///////////////////
class SEMTrajectory : public G4VTrajectory
///////////////////
{

//--------
   public:
//--------

// Constructor/Destrcutor

   SEMTrajectory();

   SEMTrajectory(const G4Track* aTrack);
   SEMTrajectory(SEMTrajectory &);
   virtual ~SEMTrajectory();

// Operators
   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const SEMTrajectory& right) const
   {return (this==&right);}

// Get/Set functions
   inline G4int GetTrackID() const
   { return fTrackID; }
   inline G4int GetParentID() const
   { return fParentID; }
   inline G4String GetParticleName() const
   { return ParticleName; }
   inline G4double GetCharge() const
   { return PDGCharge; }
   inline G4int GetPDGEncoding() const
   { return PDGEncoding; }
   inline G4ThreeVector GetInitialMomentum() const
   { return InitialMomentum; }

// Other member functions
   virtual void ShowTrajectory() const;
   virtual void ShowTrajectory(std::ostream& o) const;
   virtual void DrawTrajectory(G4int i_mode=0) const;
   virtual void AppendStep(const G4Step* aStep);
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

//---------
   private:
//---------

  SEMTrajectoryPointContainer* positionRecord;
  G4int fTrackID;
  G4int fParentID;
  G4ParticleDefinition* fpParticleDefinition;
  G4String ParticleName;
  G4double PDGCharge;
  G4int    PDGEncoding;
// FIXME not initialized !!!
  G4ThreeVector InitialMomentum;
  G4int	logenergy;
  G4double	startenergy;
  G4double  newpolyenergy;
  
//---------
   public:
//---------
   virtual int GetPointEntries() const
   { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const
   { return (*positionRecord)[i]; }
};

extern G4Allocator<SEMTrajectory> myTrajectoryAllocator;

inline void* SEMTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void SEMTrajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((SEMTrajectory*)aTrajectory);
}

#endif
