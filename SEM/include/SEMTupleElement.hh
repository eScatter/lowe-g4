#ifndef SEMTupleElement_h
#define SEMTupleElement_h

#include "globals.hh"
#include "G4ThreeVector.hh"

class SEMTupleElement
{
  public:
    G4double        toten;   // Total energy q*V+U_kin
    G4double        kinen;   // Kinetic energy U_kin
    G4double        depositen;   // Deposited energy in the step corresponding to this hit
    G4double        sort;                 // This is the sort key; copy one of the items (must be doubles) to this entry to sort on it
    G4ThreeVector   momentumdirection;    // Momentum direction in the global coordinate system
    G4ThreeVector   localposition;        // Position in the local coordinate system of the detector
    G4ThreeVector   worldposition;        // Position in the global coordinate system
    G4ThreeVector   vertexposition;       // Position where the particle was created (not necessarily the gun!)
    string          lv;                   // Logical volume where the particle was created
    string          process;              // Stores the process that created the particle
    string          category;             // Stores the category of the particle
    G4int           TrackID;
    G4int           ParentID;
    G4int           EventID;
    G4double        Time;
    G4double        Depth,Rho;

    inline bool operator<(const SEMTupleElement &rhs) const {return (this->sort < rhs.sort);}
};

#endif //SEMTupleElement_h
