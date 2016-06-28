#ifndef SEMBulkDet_h
#define SEMBulkDet_h 1

#include "SEMGeneralDetector.hh"
#include "SEMBulkDetHit.hh"
#include "SEMTuple.hh"
#include "SEMTupleElement.hh"
#include "SEMCountsPerCategory.hh"

class SEMBulkDetMessenger;

class SEMBulkDet : public SEMGeneralDetector
{

  // This is a detector that is coupled to a "Plane" type detector.

  // The energy that gets stored in the "kinen" variable of every hit comes either from the coupled Plane detector if it exists.
  // If no such particle is found, the particle is characterized as SE3, provided that during the event at least one particle 
  // passed throught the Plane detector. The stored energy is then the energy at the moment of impact on the detector.
  // If no particle passed through the coupled Plane detector, while there is still a hit on this detector a warning is issued
  // and the current energy is simply stored. 

  public:
      SEMBulkDet(G4String name);
      virtual ~SEMBulkDet();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      SEMBulkDetMessenger           *messenger;
      SEMBulkDetHitsCollection      *hitsCollection;
      G4int                     HCID;
      SEMTuple<SEMTupleElement> DetTuple;
      G4bool                    detectAtEnd,multipleHits,detectTotal;
      G4int                     fPE;
      G4double                  energyCut;

  public:
      inline void ClearTuple() {DetTuple.Clear();}
      inline SEMTuple<SEMTupleElement> *GetTuple() {return &DetTuple;}

      inline void SetMultipleHits(G4bool val) { multipleHits = val; };
      inline G4bool GetMultipleHits() const { return multipleHits; };

      inline void SetDetectAtEnd(G4bool val) { detectAtEnd = val; };
      inline G4bool GetDetectAtEnd() const { return detectAtEnd; };

      inline void SetDetectTotal(G4bool val) { detectTotal = val; };
      inline G4bool GetDetectTotal() const { return detectTotal; };

      inline void SetEnergyCut(G4double val) { energyCut = val; };
      inline G4bool GetEnergyCut() const { return energyCut; };

      void OutputCounts(string filename);
      void OutputHits(string filename);
      void OutputDeposits(string filename);
      void OutputEnergyHistogram(string filename,G4double binsize);
      void OutputEnergyPDF(string filename);
      void OutputRhoPDF(string filename);
      void OutputAngleHistogram(string filename,G4ThreeVector normal,G4int nbin);
      void OutputTimeHistogram(string filename,G4double binsize);

      inline void SetPE(G4int PE) { fPE = PE; }
      inline G4int GetPE() { return fPE; }

      inline SEMCountsPerCategory GetCounts() {return GetDetectorCounts(DetTuple);};
};

#endif

