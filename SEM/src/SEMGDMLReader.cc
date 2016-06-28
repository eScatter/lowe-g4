#include "SEMGDMLReader.hh"

void SEMGDMLReader::ExtensionRead(const xercesc::DOMElement* const extElement)
{
   G4cout << "G4GDML: Reading GDML extension..." << G4endl;

   for (xercesc::DOMNode* iter = extElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="detectorlist") { DetectorListRead(child); } else
      {
        G4String error_msg = "Unknown tag in structure: " + tag;
        G4Exception("SEMGDMLReader::ExtensionRead()",
                    "ReadError", FatalException, error_msg);
      }
   }
}

void SEMGDMLReader::VolumeRead(const xercesc::DOMElement* const volumeElement)
{
   G4cerr << "Calling SEMGDMLReader::VolumeRead" << G4endl;
   G4VSolid* solidPtr = 0;
   G4Material* materialPtr = 0;
   G4VisAttributes* attrPtr = 0;
   G4GDMLAuxListType auxList;

   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   const G4String name = Transcode(volumeElement->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = volumeElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="auxiliary")
        { auxList.push_back(AuxiliaryRead(child)); } else
      if (tag=="materialref")
        { materialPtr = GetMaterial(GenerateName(RefRead(child),true)); } else
      if (tag=="solidref")
        { solidPtr = GetSolid(GenerateName(RefRead(child))); } else
      if (tag == "color") 
        { attrPtr = ColorRead(child); } 
   }

   pMotherLogical = new G4LogicalVolume(solidPtr,materialPtr,
                                        GenerateName(name),0,0,0);
   pMotherLogical->SetVisAttributes(attrPtr);

   if (!auxList.empty()) { auxMap[pMotherLogical] = auxList; }

   Volume_contentRead(volumeElement);
   G4cerr << "Exiting SEMGDMLReader::VolumeRead" << G4endl;
}

G4VisAttributes* SEMGDMLReader::ColorRead(const xercesc::DOMElement* const colorElement)
{
   G4String name;
   G4VisAttributes* color = 0;
   G4double r=0.0, g=0.0, b=0.0, a=0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = colorElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
        { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="R")
        { r = eval.Evaluate(attValue); } else
      if (attName=="G")
        { g = eval.Evaluate(attValue); } else
      if (attName=="B")
        { b = eval.Evaluate(attValue); } else
      if (attName=="A")
        { a = eval.Evaluate(attValue); }
   }

   G4cout << "Color attribute (R,G,B,A) is: "
          << r << ", " << g << ", " << b << ", " << a << " !" << G4endl;
   color = new G4VisAttributes(G4Color(r,g,b,a));
   return color;
}

void SEMGDMLReader::DetectorListRead(const xercesc::DOMElement* const detectorListElement)
{
   G4cout << "SEMGDML: Reading detectors ..." << G4endl;

   for (xercesc::DOMNode* iter = detectorListElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="detector") { DetectorRead(child); }
      else
      {
        G4String error_msg = "Unknown tag in structure: " + tag;
        G4Exception("SEMGDMLReader::FieldListRead()",
                    "ReadError", FatalException, error_msg);
      }
   }
}

void SEMGDMLReader::DetectorRead(const xercesc::DOMElement* const element)
{
   G4cout << "SEMGDML: Reading detector..." << G4endl;

   G4String name,type,volume,coupledname;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
        { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
      const G4String attName  = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name")  { name = attValue; } else
      if (attName=="volume"){ volume = attValue; } else
      if (attName=="type")  { type = attValue; }
   }

   // Store the data for the detectors .... TBD
   G4int ndet = detectorlist.size();
   G4cerr << "Number of stored detectors = " << ndet << G4endl;
   CADDetectorData* det = new CADDetectorData;
   
   det->SetName(name);
   det->SetVolume(volume);
   det->SetCoupledName(coupledname);
   CADDtype detectortype=fNone;
   if( type != "" ) {
      if( type=="plane" ) detectortype = fPlane; else
      if( type=="bulk" ) detectortype = fBulk; else
      if( type=="counter" ) detectortype = fCounter; else
      if( type=="user0" ) detectortype = fUser0; else
      if( type=="user1" ) detectortype = fUser1; else
      if( type=="user2" ) detectortype = fUser2; else
      if( type=="user3" ) detectortype = fUser3; else
      if( type=="user4" ) detectortype = fUser4;
    }
    det->SetType(detectortype);
    detectorlist.push_back(det);

}

