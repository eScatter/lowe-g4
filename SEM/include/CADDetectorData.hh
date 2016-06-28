//
#ifndef CADDETECTORDATA_H
#define CADDETECTORDATA_H 1

// G4 includes
#include "globals.hh"

// Definitions for detectortypes
enum CADDtype {fNone = 0,fPlane = 10,fBulk = 20,fCounter = 30, fUser0 = 100,fUser1 = 110, fUser2 = 120,fUser3 = 130,fUser4 = 140};

class CADDetectorData
{
public:
   inline void SetName(G4String n) {name = n;}
   inline G4String GetName() {return name;}

   inline void SetVolume(G4String vn) {logvolume = vn;}
   inline G4String GetVolume() {return logvolume;}

   inline void SetCoupledName(G4String cn) {coupledname = cn;}
   inline G4String GetCoupledName() {return coupledname;}
    
   inline void SetType(CADDtype t) {type = t;}
   inline CADDtype GetType() {return type;}

private:
   G4String name;
   G4String logvolume;
   G4String coupledname;
   CADDtype type;
};

#endif // CADDETECTORDATA_H

