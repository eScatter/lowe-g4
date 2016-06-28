#ifndef SEMDetectorConstruction_h
#define SEMDetectorConstruction_h 1

#include "globals.hh"
#include <vector>
#include "G4VUserDetectorConstruction.hh"

//EB: needed for GDML input:
#include "G4GDMLParser.hh" 
#include "SEMGDMLReader.hh"
//EB
#include "CADDetectorData.hh"

class SEMDetectorConstructionMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

using namespace std;

class SEMDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     SEMDetectorConstruction();
     SEMDetectorConstruction(const G4String& filename);
     virtual ~SEMDetectorConstruction();

  public:
  
     virtual G4VPhysicalVolume* Construct();
     G4LogicalVolume* FindLogicalVolume(const G4String&);
     void DumpGDMLGeometry();
     
  private:

     void DumpGeometricalTree(G4VPhysicalVolume* aVolume,ostream& out,G4int depth=0);

     G4GDMLParser *parser;
     SEMGDMLReader *reader;

     G4VPhysicalVolume      *fWorld;       // pointer to the physical envelope
     
     SEMDetectorConstructionMessenger *messenger;  // pointer to the messenger
     
     vector < CADDetectorData* >  *detectorlist;

  public:
     inline G4int GetNumberOfDetectors() { return detectorlist->size(); };
     inline G4String GetDetectorName(G4int i) { return ((*detectorlist)[i])->GetName(); };
     inline G4String GetDetectorCoupledName(G4int i) { return ((*detectorlist)[i])->GetCoupledName(); };
     inline G4int GetDetectorType(G4int i) { return ((*detectorlist)[i])->GetType(); };

     void setMaterial (G4String volumename, G4String materialname);
     void setFraction (G4String materialname, G4int index, G4double fraction);
     void setGasDensity (G4String materialname, G4double density);
     G4double getGasDensity (G4String materialname);
     void PrintTree(ostream& out);
};
#endif
