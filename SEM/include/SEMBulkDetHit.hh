#ifndef SEMBulkDetHit_h
#define SEMBulkDetHit_h 1

#include "globals.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class SEMBulkDetHit : public G4VHit
{
  public:

      SEMBulkDetHit();
      virtual ~SEMBulkDetHit();
      SEMBulkDetHit(const SEMBulkDetHit &right);
      const SEMBulkDetHit& operator=(const SEMBulkDetHit &right);
      int operator==(const SEMBulkDetHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      inline float x();
      inline float y();

      virtual void Draw();
      virtual void Print();

  private:
      G4double time;
      G4ThreeVector localPos;
      G4ThreeVector worldPos;
      G4ThreeVector momentumdirection;
      G4ThreeVector vertexPos;
      G4double      kinenergy;
      G4double      totalenergy;
      G4double      depositedenergy;
      G4String      lv;
      G4String      process;
      G4String      category;
      G4int         ParentID;
      G4int         TrackID;
      G4int         EventID;
      G4double      Depth,Rho;

  public:
      inline void SetTime(G4double t) { time = t; }
      inline G4double GetTime() const { return time; }
      inline void SetLocalPos(G4ThreeVector xyz) { localPos = xyz; }
      inline G4ThreeVector GetLocalPos() const { return localPos; }
      inline void SetVertexPos(G4ThreeVector xyz) { vertexPos = xyz; }
      inline G4ThreeVector GetVertexPos() const { return vertexPos; }
      inline void SetWorldPos(G4ThreeVector xyz) { worldPos = xyz; }
      inline G4ThreeVector GetWorldPos() const { return worldPos; }
      inline void SetMomentumDirection(G4ThreeVector xyz) { momentumdirection = xyz; }
      inline G4ThreeVector GetMomentumDirection() const { return momentumdirection; }
      inline void SetTotalEnergy(G4double E) { totalenergy = E; }
      inline G4double GetTotalEnergy() const { return totalenergy; }
      inline void SetKinEnergy(G4double kin) { kinenergy = kin; }
      inline G4double GetKinEnergy() const { return kinenergy; }
      inline void SetDepositedEnergy(G4double E) { depositedenergy = E; }
      inline G4double GetDepositedEnergy() const { return depositedenergy; }
      inline void SetLVAtVertex(const G4LogicalVolume* lvatvertex) {lv = lvatvertex->GetName();}
      inline G4String GetLVAtVertex() { return lv; }
      inline void SetCreatorProcess(const G4String name) {process=name; }
      inline G4String GetCreatorProcess() {return process; }
      inline void SetCategory(const G4String name) {category=name; }
      inline G4String GetCategory() {return category; }
      inline G4int GetParentID() {return ParentID;}
      inline void SetParentID(G4int ID) {ParentID=ID;}
      inline G4int GetTrackID() {return TrackID;}
      inline void SetTrackID(G4int ID) {TrackID=ID;}
      inline G4int GetEventID() {return EventID;}
      inline void SetEventID(G4int ID) {EventID=ID;}
      inline void SetDepth(G4double d) { Depth = d; }
      inline G4double GetDepth() const { return Depth; }
      inline void SetRho(G4double d) { Rho = d; }
      inline G4double GetRho() const { return Rho; }

};

typedef G4THitsCollection<SEMBulkDetHit> SEMBulkDetHitsCollection;

extern G4Allocator<SEMBulkDetHit> SEMBulkDetHitAllocator;

inline void* SEMBulkDetHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)SEMBulkDetHitAllocator.MallocSingle();
  return aHit;
}

inline void SEMBulkDetHit::operator delete(void* aHit)
{
  SEMBulkDetHitAllocator.FreeSingle((SEMBulkDetHit*) aHit);
}

#endif


