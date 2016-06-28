#ifndef SEMGDMLReader_H
#define SEMGDMLReader_H 1

#include <map>
#include "G4GDMLReadStructure.hh"
#include "CADDetectorData.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

class SEMGDMLReader : public G4GDMLReadStructure
{

 public:

   void ExtensionRead(const xercesc::DOMElement* const element);
   void VolumeRead(const xercesc::DOMElement* const element);
  
   inline std::vector<CADDetectorData*>* GetDetectorList() {return &detectorlist;}

 protected:

   void DetectorListRead(const xercesc::DOMElement* const element);
   void ElectrodeRead(const xercesc::DOMElement* const element,G4int n);
   void DetectorRead(const xercesc::DOMElement* const element);
   G4VisAttributes* ColorRead(const xercesc::DOMElement* const element);

 private:

   std::vector<CADDetectorData*> detectorlist;
};

#endif
