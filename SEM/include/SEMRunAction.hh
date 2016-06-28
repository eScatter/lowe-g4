#ifndef SEMRunAction_h
#define SEMRunAction_h 1

#include "G4UserRunAction.hh"

#include <vector>
#include "SEMTuple.hh"
#include "SEMTupleElement.hh"
#include "G4Colour.hh"

using namespace std;

class G4Run;
class SEMScanMessenger;
class SEMCountsPerCategory;

class SEMRunAction : public G4UserRunAction
{
  public:
    SEMRunAction();
   ~SEMRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
 
  private:
    G4int PE;
    void CountCategories(ofstream *output,SEMTuple<SEMTupleElement> tuple);

  public:
    virtual void MapOn(G4String,G4int,G4int);
    virtual void ScanOn(G4int);
    virtual void ScanOn(G4int,G4String);
    virtual void ImageOn(G4int,G4String);

  private: 
    SEMScanMessenger* scanMessenger;
    G4int scanNx;
    G4int scanNy;
    G4double scanLx;
    G4double scanLy;
    G4ThreeVector scanXdir;
    G4ThreeVector scanYdir;
    G4ThreeVector scanOrigin;
    G4bool        validscan,validimage;
    G4bool        maptrackend;
    G4Colour      mapcolour;
    G4String      scanPrefix;
    G4ThreeVector mapRotation;

    vector < G4ThreeVector > fPosLineScan;
    vector < vector < SEMCountsPerCategory > > fLineScan;
    vector < vector < vector < SEMCountsPerCategory > > > fImage;

    void WritePgm(G4String,vector<vector<SEMCountsPerCategory > >);
    void WritePgm(G4String,vector<vector<SEMCountsPerCategory > >,G4int,G4bool,G4bool,G4double*, G4double*);

  public:
    inline void SetScanNx(G4int num) { scanNx = num; validscan = false; validimage = false; }
    inline G4int GetScanNx() const { return scanNx; }
    inline void SetScanNy(G4int num) { scanNy = num; validimage = false; }
    inline G4int GetScanNy() const { return scanNy; }
    inline void SetScanLx(G4double num) { scanLx = num; }
    inline G4double GetScanLx() const { return scanLx; }
    inline void SetScanLy(G4double num) { scanLy = num; }
    inline G4double GetScanLy() const { return scanLy; }
    inline void SetScanXdir(G4ThreeVector vec) { scanXdir = vec/vec.mag(); }
    inline G4ThreeVector GetScanXdir() const { return scanXdir; }
    inline void SetScanYdir(G4ThreeVector vec) { scanYdir = vec/vec.mag(); }
    inline G4ThreeVector GetScanYdir() const { return scanYdir; }
    inline void SetScanOrigin(G4ThreeVector vec) { scanOrigin = vec; }
    inline G4ThreeVector GetScanOrigin() const { return scanOrigin; }
    inline void SetScanPrefix(G4String prefix) { scanPrefix = prefix; }
    inline G4String GetScanPrefix() const { return scanPrefix; }
    inline void SetMapRotation(G4ThreeVector vec) {mapRotation = vec;}

    void OutputLinescan(G4String filename,G4String detectorname);
    void OutputPgm(G4String filename,G4String detectorname,G4int category,G4bool autoscale, G4double min, G4double max);
    void OutputASCII(G4String filename,G4String detectorname);

    inline G4bool GetMapTrackEnd() {return maptrackend;}
    inline G4Colour GetMapColour() {return mapcolour;}
    inline void SetMapColour(G4Colour mc) {mapcolour = mc;}

};

#endif
