// Oct 16 2008 Created

#ifndef SEMCounter_h
#define SEMCounter_h 1

#include "SEMGeneralDetector.hh"
#include "SEMCounterHit.hh"
#include "SEMTuple.hh"
#include "SEMTupleElement.hh"
#include "G4ThreeVector.hh"
#include "SEMCountsPerCategory.hh"

class SEMCounterMessenger;

class SEMCounter : public SEMGeneralDetector
{

  // This detector stores information about the total deposited energy for each event

  public:
      SEMCounter(G4String name);
      virtual ~SEMCounter();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      SEMCounterMessenger*        messenger;
      SEMCounterHitsCollection*   hitsCollection;
      G4int                     HCID;
      G4bool                    directional;
      G4bool                    AllowMultipleHits;
      G4bool                    CountEHPairs;
      G4ThreeVector             direction;
      SEMTuple<SEMTupleElement> CounterTuple;
      G4int                     fPE,fSE,fBSE,nhit;
      G4double                  EHit;
      G4int                     lastTID;
      G4int                     verboseLevel;
      G4double                  DepositedTillNow;


  public:
      inline void ClearTuple() {CounterTuple.Clear(); fSE=fBSE=0; }
      inline SEMTuple<SEMTupleElement> *GetTuple() {return &CounterTuple;}

      void OutputCounts(string filename);
      void OutputHits(string filename);
      void OutputEnergyHistogram(string filename,G4double binsize);
      void OutputEnergyPDF(string filename);

      inline void SetPE(G4int PE) { fPE = PE; }
      inline G4int GetPE() { return fPE; }
      inline G4int GetSE() { return fSE; };
      inline G4int GetBSE() { return fBSE; };

      inline SEMCountsPerCategory GetCounts() {return GetDetectorCounts(CounterTuple);};
      inline void SetVerbose(int v) {verboseLevel=v;}
};

#endif

