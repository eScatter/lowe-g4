// May 21 2007: Added new version of angular histogram using energy selection

#ifndef SEMPlane_h
#define SEMPlane_h 1

#include "SEMGeneralDetector.hh"
#include "SEMPlaneHit.hh"
#include "SEMTuple.hh"
#include "SEMTupleElement.hh"
#include "G4ThreeVector.hh"
#include "SEMCountsPerCategory.hh"

class SEMPlaneMessenger;

class SEMPlane : public SEMGeneralDetector
{

  // This detector stores information about the FIRST impact by any given particle. 

  public:
      SEMPlane(G4String name);
      virtual ~SEMPlane();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      SEMPlaneMessenger*        messenger;
      SEMPlaneHitsCollection*   hitsCollection;
      G4int                     HCID;
      G4bool                    directional;
      G4bool                    AllowMultipleHits;
      G4bool                    CountEHPairs;
      G4ThreeVector             direction;
      SEMTuple<SEMTupleElement> PlaneTuple;
      G4int                     fPE,fSE,fBSE;


  public:
      inline void SetDirectional(G4bool val) { directional = val; };
      inline G4bool GetDirectional() const { return directional; };
      inline void SetDirection(G4ThreeVector vec) { direction = vec; };
      inline void SetDirection(G4double Theta, G4double Phi) {G4cout << "Direction: theta = " << Theta/deg << " deg   phi = " << Phi/deg << " deg" << G4endl;  direction = G4ThreeVector(sin(Theta)*cos(Phi),sin(Theta)*sin(Phi),cos(Theta)); };

      inline void SetMultipleHits(G4bool val) { AllowMultipleHits = val; };
      inline G4bool GetMultipleHits() const { return AllowMultipleHits; };

      inline void SetCountEHPairs(G4bool val) { CountEHPairs = val; };
      inline G4bool GetCountEHPairs() const { return CountEHPairs; };

      inline void ClearTuple() {PlaneTuple.Clear(); fSE=fBSE=0; }
      inline SEMTuple<SEMTupleElement> *GetTuple() {return &PlaneTuple;}

      void OutputCounts(string filename);
      void OutputHits(string filename);
      void OutputDeposits(string filename,G4double mindeposit);
      void OutputEnergyHistogram(string filename,G4double binsize);
      void OutputEnergyPDF(string filename);
      void OutputGammaPDF(string filename);
      void OutputRhoHistogram(string filename,G4double binsize);
      void OutputRhoPDF(string filename);
      void OutputAngleHistogram(string filename,G4ThreeVector normal,G4int nbin);

      inline void SetPE(G4int PE) { fPE = PE; }
      inline G4int GetPE() { return fPE; }
      inline G4int GetSE() { return fSE; };
      inline G4int GetBSE() { return fBSE; };

      inline SEMCountsPerCategory GetCounts() {return GetDetectorCounts(PlaneTuple);};
};

#endif

