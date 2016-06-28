//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4GDMLReadMaterials.cc 70764 2013-06-05 12:54:37Z gcosmo $
// GEANT4 tag $ Name:$
//
// class G4GDMLReadMaterials Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLReadMaterials.hh"

#include <vector>
#include <iostream>
using namespace std;
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

G4GDMLReadMaterials::G4GDMLReadMaterials() : G4GDMLReadDefine()
{
}

G4GDMLReadMaterials::~G4GDMLReadMaterials()
{
}

G4double
G4GDMLReadMaterials::AtomRead(const xercesc::DOMElement* const atomElement)
{
   G4double value = 0.0;
   G4double unit = g/mole;

   const xercesc::DOMNamedNodeMap* const attributes
         = atomElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::AtomRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return value;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") { value = eval.Evaluate(attValue); } else
      if (attName=="unit")  { unit = eval.Evaluate(attValue); }
   }

   return value*unit;
}

G4int G4GDMLReadMaterials::
CompositeRead(const xercesc::DOMElement* const compositeElement,G4String& ref)
{
   G4int n = 0;

   const xercesc::DOMNamedNodeMap* const attributes
         = compositeElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::CompositeRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return n;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="n")  { n = eval.EvaluateInteger(attValue); } else
      if (attName=="ref") { ref = attValue; }
   }

   return n;
}

G4double G4GDMLReadMaterials::DRead(const xercesc::DOMElement* const DElement)
{
   G4double value = 0.0;
   G4double unit = g/cm3;

   const xercesc::DOMNamedNodeMap* const attributes
         = DElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::DRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return value;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") { value = eval.Evaluate(attValue); } else
      if (attName=="unit")  { unit = eval.Evaluate(attValue); }
   }

   return value*unit;
}

G4double G4GDMLReadMaterials::PRead(const xercesc::DOMElement* const PElement)
{
   G4double value = STP_Pressure;
   G4double unit = hep_pascal;

   const xercesc::DOMNamedNodeMap* const attributes = PElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::PRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return value;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") { value = eval.Evaluate(attValue); } else
      if (attName=="unit")  { unit = eval.Evaluate(attValue); }
   }

   return value*unit;
}

G4double G4GDMLReadMaterials::TRead(const xercesc::DOMElement* const TElement)
{
   G4double value = STP_Temperature;
   G4double unit = kelvin;

   const xercesc::DOMNamedNodeMap* const attributes = TElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::TRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return value;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") { value = eval.Evaluate(attValue); } else
      if (attName=="unit")  { unit = eval.Evaluate(attValue); }
   }

   return value*unit;
}

G4double G4GDMLReadMaterials::MEERead(const xercesc::DOMElement* const PElement)
{
   G4double value = -1;
   G4double unit = eV;

   const xercesc::DOMNamedNodeMap* const attributes = PElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::MEERead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return value;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") { value = eval.Evaluate(attValue); } else
      if (attName=="unit")  { unit = eval.Evaluate(attValue); }
   }

   return value*unit;
}

void G4GDMLReadMaterials::
ElementRead(const xercesc::DOMElement* const elementElement) 
{
   G4String name;
   G4String formula;
   G4double a = 0.0;
   G4double Z = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = elementElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::ElementRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); } else
      if (attName=="formula") { formula = attValue; } else
      if (attName=="Z") { Z = eval.Evaluate(attValue); }
   }

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = elementElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::ElementRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="atom") { a = AtomRead(child); }  else
      if (tag=="fraction") { nComponents++; }
   }

   if (nComponents>0)
   {
     MixtureRead(elementElement,
                 new G4Element(Strip(name),formula,nComponents));
   }
   else
   {
     new G4Element(Strip(name),formula,Z,a);
   }
}

G4double G4GDMLReadMaterials::
FractionRead(const xercesc::DOMElement* const fractionElement, G4String& ref)
{
   G4double n = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = fractionElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::FractionRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return n;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="n")   { n = eval.Evaluate(attValue); } else
      if (attName=="ref") { ref = attValue; }
   }

   return n;
}

void G4GDMLReadMaterials::
IsotopeRead(const xercesc::DOMElement* const isotopeElement)
{
   G4String name;
   G4int Z = 0;
   G4int N = 0;
   G4double a = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = isotopeElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount;attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::IsotopeRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); } else
      if (attName=="Z") { Z = eval.EvaluateInteger(attValue); } else
      if (attName=="N") { N = eval.EvaluateInteger(attValue); }
   }

   for (xercesc::DOMNode* iter = isotopeElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::IsotopeRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="atom")  { a = AtomRead(child); }
   }

   new G4Isotope(Strip(name),Z,N,a);
}

void G4GDMLReadMaterials::
MaterialRead(const xercesc::DOMElement* const materialElement)
{
   G4String name;
   G4double Z = 0.0;
   G4double a = 0.0;
   G4double D = 0.0;
   G4State state = kStateUndefined;
   G4double T = STP_Temperature;
   G4double P = STP_Pressure;
   //EB: Added
   G4double   mwf = 0.0;
   G4double   fermi = 0.0;
   G4double   bb = 0.0;
   G4double   deltaphi = 0.0;
   G4double   affinity = 0.0;
   G4double   bandgap = 0.0;
   G4int      conductortype = 0;
   G4double   resistivity = 0.0;
   G4double   soundvelocity = 0.0;
   G4double   defpotential = 0.0;
   G4double   lattice = 0.0;
   G4String   formula = "";
   vector < double > gasfractions;
   G4double   gasdensity = 0.0;
   G4bool hasproperties = false;
   // EB: End of addition
   G4double MEE = -1.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = materialElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::MaterialRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); } else
      if (attName=="Z") { Z = eval.Evaluate(attValue); } else
      if (attName=="state")
      {
         if (attValue=="solid")  { state = kStateSolid;  } else
         if (attValue=="liquid") { state = kStateLiquid; } else
         if (attValue=="gas")    { state = kStateGas; }
      } else
      // EB: Added properties
      if (attName=="formula") { formula = attValue;} else
      if (attName=="workfunction") { mwf = eval.Evaluate(attValue)*eV; if (0!=mwf) hasproperties=true;} else
      if (attName=="fermienergy") { fermi = eval.Evaluate(attValue)*eV; if (0!=fermi) hasproperties=true;} else
      if (attName=="bandbending") { bb = eval.Evaluate(attValue)*eV; if (0!=bb) hasproperties=true;} else
      if (attName=="deltaphi") { deltaphi = eval.Evaluate(attValue)*eV; if (0!=deltaphi) hasproperties=true;} else
      if (attName=="affinity") { affinity = eval.Evaluate(attValue)*eV; if (0!=affinity) hasproperties=true;} else
      if (attName=="bandgap") { bandgap = eval.Evaluate(attValue)*eV; if (0!=bandgap) hasproperties=true;} else
      if (attName=="conductortype")
      {
         if (attValue=="insulator")  { conductortype = 2;  } else
         if (attValue=="semi") { conductortype = 1; } 
         if (0!=conductortype) hasproperties=true;
      } else
      if (attName=="resistivity") { resistivity = eval.Evaluate(attValue); if (0!=resistivity) hasproperties=true;} else
      if (attName=="soundvelocity") { soundvelocity = eval.Evaluate(attValue); if (0!=soundvelocity) hasproperties=true;} else
      if (attName=="defpotential") { defpotential = eval.Evaluate(attValue)*eV; if (0!=defpotential) hasproperties=true;} else
      if (attName=="lattice") { lattice = eval.Evaluate(attValue); if (0!=lattice) hasproperties=true;} else
      if (attName=="gasfractions")
      {
         istringstream is;
         is.str(attValue);
         double fraction;
         while (!is.eof()) {
            is >> fraction;
            gasfractions.push_back(fraction);
         }
         hasproperties=true;
      } else
      if (attName=="gasdensity") { gasdensity = eval.Evaluate(attValue); if (0!=gasdensity) hasproperties=true;} 
      // EB: End of additional properties
   }

   size_t nComponents = 0;

   for (xercesc::DOMNode* iter = materialElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::MaterialRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="atom") { a = AtomRead(child); } else
      if (tag=="Dref") { D = GetQuantity(GenerateName(RefRead(child))); } else
      if (tag=="Pref") { P = GetQuantity(GenerateName(RefRead(child))); } else
      if (tag=="Tref") { T = GetQuantity(GenerateName(RefRead(child))); } else
      if (tag=="MEEref") { MEE = GetQuantity(GenerateName(RefRead(child))); } else
      if (tag=="D") { D = DRead(child); } else
      if (tag=="P") { P = PRead(child); } else
      if (tag=="T") { T = TRead(child); } else
      if (tag=="MEE") { MEE = MEERead(child); } else
      if (tag=="fraction" || tag=="composite")  { nComponents++; }
   }

   G4Material* material =  0;

   if (nComponents==0)
   {
     material = new G4Material(Strip(name),Z,a,D,state,T,P);
   }
   else
   {
     material = new G4Material(Strip(name),D,nComponents,state,T,P);
     MixtureRead(materialElement, material);
   }

   // EB: Added 
   material->SetChemicalFormula(formula);
   if (hasproperties) {
      G4MaterialPropertiesTable* matprop = material->GetMaterialPropertiesTable();
      if (!matprop)
      {
        material->SetMaterialPropertiesTable(
                  matprop = new G4MaterialPropertiesTable());
      }
      if (mwf!=0) matprop->AddConstProperty("WORKFUNCTION",mwf);
      if (fermi!=0) matprop->AddConstProperty("FERMIENERGY",fermi);
      if (bb!=0) matprop->AddConstProperty("BANDBENDING",bb);
      if (deltaphi!=0) matprop->AddConstProperty("DELTAPHI",deltaphi);
      if (affinity!=0) matprop->AddConstProperty("AFFINITY",affinity);
      if (bandgap!=0) matprop->AddConstProperty("BANDGAP",bandgap);
      if (conductortype!=0) matprop->AddConstProperty("CONDUCTORTYPE",conductortype);
      if (resistivity!=0) matprop->AddConstProperty("RESISTIVITY",resistivity);
      if (soundvelocity!=0) matprop->AddConstProperty("SOUNDVELOCITY",soundvelocity);
      if (defpotential!=0) matprop->AddConstProperty("DEFPOTENTIAL",defpotential);
      if (lattice!=0) matprop->AddConstProperty("LATTICE",lattice);
      char property[10]="";
      for (size_t i=0;i<gasfractions.size();i++) {
         sprintf (property,"FRACTION%d",(int) i);
         matprop->AddConstProperty(property,gasfractions[i]);
      }
      if (gasdensity!=0) matprop->AddConstProperty("GASDENSITY",gasdensity);
   }
   // EB: End of addition
   if (MEE != -1)  // ionisation potential (mean excitation energy)
   {
     material->GetIonisation()->SetMeanExcitationEnergy(MEE);
   }

   for (xercesc::DOMNode* iter = materialElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::MaterialRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="property") { PropertyRead(child,material); }
   }
}

void G4GDMLReadMaterials::
MixtureRead(const xercesc::DOMElement *const mixtureElement, G4Element *element)
{
   for (xercesc::DOMNode* iter = mixtureElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::MixtureRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="fraction")
      {
         G4String ref;
         G4double n = FractionRead(child,ref);
         element->AddIsotope(GetIsotope(GenerateName(ref,true)),n);
      }
   }
}

void G4GDMLReadMaterials::
MixtureRead(const xercesc::DOMElement *const mixtureElement,
            G4Material *material)
{
   for (xercesc::DOMNode* iter = mixtureElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::MixtureRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="fraction")
      {
         G4String ref;
         G4double n = FractionRead(child,ref);
         
         G4Material *materialPtr = GetMaterial(GenerateName(ref,true), false);
         G4Element *elementPtr = GetElement(GenerateName(ref,true), false);

         if (materialPtr != 0) { material->AddMaterial(materialPtr,n); } else
         if (elementPtr != 0)  { material->AddElement(elementPtr,n); }

         if ((materialPtr == 0) && (elementPtr == 0))
         {
            G4String error_msg = "Referenced material/element '"
                               + GenerateName(ref,true) + "' was not found!";
            G4Exception("G4GDMLReadMaterials::MixtureRead()", "InvalidSetup",
                        FatalException, error_msg);   
         }
      } 
      else if (tag=="composite")
      {
         G4String ref;
         G4int n = CompositeRead(child,ref);

         G4Element *elementPtr = GetElement(GenerateName(ref,true));
         if (elementPtr == 0)
         {
            G4String error_msg = "Referenced material/element '"
                               + GenerateName(ref,true) + "' was not found!";
            G4Exception("G4GDMLReadMaterials::MixtureRead()", "InvalidSetup",
                        FatalException, error_msg);   
         }
         material->AddElement(elementPtr,n);
      }
   }
}

void G4GDMLReadMaterials::
PropertyRead(const xercesc::DOMElement* const propertyElement,
             G4Material* material)
{
   G4String name;
   G4String ref;
   G4GDMLMatrix matrix;

   const xercesc::DOMNamedNodeMap* const attributes
         = propertyElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLReadMaterials::PropertyRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); } else
      if (attName=="ref")  { matrix = GetMatrix(ref=attValue); }
   }

   /*
   if (matrix.GetCols() != 2)
   {
     G4String error_msg = "Referenced matrix '" + ref
            + "' should have \n two columns as a property table for material: "
            + material->GetName();
     G4Exception("G4GDMLReadMaterials::PropertyRead()", "InvalidRead",
                 FatalException, error_msg);
   }
   */

   if (matrix.GetRows() == 0) { return; }

   G4MaterialPropertiesTable* matprop=material->GetMaterialPropertiesTable();
   if (!matprop)
   {

     matprop = new G4MaterialPropertiesTable();
     material->SetMaterialPropertiesTable(matprop);
   }


   if (matrix.GetCols() == 1)  // constant property assumed
   {
     matprop->AddConstProperty(Strip(name), matrix.Get(0,0));
   }
   else  // build the material properties vector
   {
     G4MaterialPropertyVector* propvect = new G4MaterialPropertyVector();
     for (size_t i=0; i<matrix.GetRows(); i++)
     {
       propvect->InsertValues(matrix.Get(i,0),matrix.Get(i,1));
     }
     matprop->AddProperty(Strip(name),propvect);
   }
}

void G4GDMLReadMaterials::
MaterialsRead(const xercesc::DOMElement* const materialsElement)
{
   G4cout << "G4GDML: Reading materials..." << G4endl;

   for (xercesc::DOMNode* iter = materialsElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLReadMaterials::MaterialsRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      
      if (tag=="define")   { DefineRead(child);  }  else 
      if (tag=="element")  { ElementRead(child); }  else 
      if (tag=="isotope")  { IsotopeRead(child); }  else 
      if (tag=="material") { MaterialRead(child); }
      else
      {
        G4String error_msg = "Unknown tag in materials: " + tag;
        G4Exception("G4GDMLReadMaterials::MaterialsRead()", "InvalidSetup",
                    FatalException, error_msg);
      }
   }
}

G4Element* G4GDMLReadMaterials::
GetElement(const G4String& ref, G4bool verbose) const
{
   G4Element* elementPtr = G4Element::GetElement(ref,false);

   if (!elementPtr)
   {
     elementPtr = G4NistManager::Instance()->FindOrBuildElement(ref);
   }

   if (verbose && !elementPtr)
   {
     G4String error_msg = "Referenced element '" + ref + "' was not found!";
     G4Exception("G4GDMLReadMaterials::GetElement()", "InvalidRead",
                 FatalException, error_msg);
   }

   return elementPtr;
}

G4Isotope* G4GDMLReadMaterials::GetIsotope(const G4String& ref,
                                           G4bool verbose) const
{
   G4Isotope* isotopePtr = G4Isotope::GetIsotope(ref,false);

   if (verbose && !isotopePtr)
   {
     G4String error_msg = "Referenced isotope '" + ref + "' was not found!";
     G4Exception("G4GDMLReadMaterials::GetIsotope()", "InvalidRead",
                 FatalException, error_msg);
   }

   return isotopePtr;
}

G4Material* G4GDMLReadMaterials::GetMaterial(const G4String& ref,
                                             G4bool verbose) const
{
   G4Material *materialPtr = G4Material::GetMaterial(ref,false);

   if (!materialPtr)
   {
     materialPtr = G4NistManager::Instance()->FindOrBuildMaterial(ref);
   }

   if (verbose && !materialPtr)
   {
     G4String error_msg = "Referenced material '" + ref + "' was not found!";
     G4Exception("G4GDMLReadMaterials::GetMaterial()", "InvalidRead",
                 FatalException, error_msg);
   }

   return materialPtr;
}
