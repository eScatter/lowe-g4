#include "SEMDetectorConstruction.hh"
#include "SEMDetectorConstructionMessenger.hh"
#include "G4LogicalVolume.hh"

#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "SEMBulkDet.hh"
#include "SEMPlane.hh"
#include "SEMCounter.hh"
#include "G4RunManagerKernel.hh" 
#include "G4LogicalVolumeStore.hh"

#include "SEMGDMLReader.hh"

#include "CLHEP/Evaluator/Evaluator.h"

SEMDetectorConstruction::SEMDetectorConstruction()
:fWorld(0)
{
  messenger = new SEMDetectorConstructionMessenger(this);

  reader = new SEMGDMLReader;
  parser = new G4GDMLParser(reader);
  parser->Read("SEM.gdml", false); // Read the gdml input, don't bother comparing to the gdml schema as it is not yet updated
}

SEMDetectorConstruction::SEMDetectorConstruction(const G4String& filename)
:fWorld(0)
{
  messenger = new SEMDetectorConstructionMessenger(this);

  reader = new SEMGDMLReader;
  parser = new G4GDMLParser(reader);
  parser->Read(filename, false); // Read the gdml input, don't bother comparing to the gdml schema as it is not yet updated
}


SEMDetectorConstruction::~SEMDetectorConstruction()
{
  if (reader) delete reader;
  if (parser) delete parser;
}

G4VPhysicalVolume* SEMDetectorConstruction::Construct()
{

  // Retrieve the world volume from GDML processor
  fWorld = parser->GetWorldVolume();
  // Check if the world volume seems to be fine
  if ( fWorld == 0 ) {
     G4Exception("Invalid world volume, check your setup selection criteria or GDML input!","VolumeError",FatalException,"SEMDetectorConstruction(1)");
  }


 // Define the detectors if they exist
 G4LogicalVolume* lv =0;
 G4SDManager* SDman = G4SDManager::GetSDMpointer();

 //detectorlist = GDMLProcessor::GetInstance()->GetDetectorList();
 detectorlist = reader->GetDetectorList();

 G4int ndetectors = detectorlist->size();
 if (ndetectors==0) {
    G4cerr << "****** SEMDetectorConstruction: No detector definitions found, you may want to check your GDML input!" << G4endl;
 }
 // Definition of detectors
 for (int i=0; i<ndetectors; i++) {
    CADDetectorData* detector = (*detectorlist)[i];
    lv = FindLogicalVolume(detector->GetVolume());
    if (lv) {
       switch (detector->GetType())
       {
         SEMPlane *adetector;
         SEMBulkDet *bdetector;
         SEMCounter *cdetector;
         case fPlane:
            adetector = new SEMPlane(detector->GetName());
            SDman->AddNewDetector(adetector);
            lv->SetSensitiveDetector(adetector);
            adetector->SetLVname(lv->GetName());
            G4cerr << "Name of planar detector stored = " << "#" << adetector->GetName() << "#" << G4endl;
            break;
         case fCounter:
            cdetector = new SEMCounter(detector->GetName());
            SDman->AddNewDetector(cdetector);
            lv->SetSensitiveDetector(cdetector);
            cdetector->SetLVname(lv->GetName());
            G4cerr << "Name of counter stored = " << "#" << cdetector->GetName() << "#" << G4endl;
            break;
         case fBulk:
            bdetector = new SEMBulkDet(detector->GetName());
            SDman->AddNewDetector(bdetector);
            lv->SetSensitiveDetector(bdetector);
            bdetector->SetLVname(lv->GetName());
            G4cerr << "Name of detector stored = " << "#" << bdetector->GetName() << "#" << G4endl;
         G4cerr << "Material of volume " << lv->GetName() << " changed to BlackHole." << G4endl;
            break;
         default:
            G4Exception("Unknown detectortype encountered, check your GDML input!","DetectorTypeError",FatalException,"SEMDetectorConstruction(2)");
       }
    } else {
       G4String msg;
       msg = "No logical volume "+detector->GetVolume()+" found (needed for detector "+detector->GetName()+"). Check your GDML input!";
       G4Exception(msg,"DetectorTypeError",FatalException,"SEMDetectorConstruction(2)");
    }
 }

   ///////////////////////////////////////////////////////////////////////
   //
   // Retrieve Auxiliary Information (color)
   //
   G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
   std::vector<G4LogicalVolume*>::iterator lvciter;
   for( lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++ )
   { 
     G4GDMLAuxListType auxInfo = parser->GetVolumeAuxiliaryInformation(*lvciter);
     std::vector<G4GDMLAuxStructType>::const_iterator ipair = auxInfo.begin();
     for( ipair = auxInfo.begin(); ipair != auxInfo.end(); ipair++ )
     {
       HepTool::Evaluator eval;
       G4String str=ipair->type;
       G4double value=eval.evaluate(ipair->value);
       G4double A,R,G,B;
       G4double tmp;
       G4VisAttributes* color;
       tmp = value - floor(value);
       R = tmp*2.0;
       value = (value-R*0.5)/1000.0;
       tmp = value - floor(value);
       G = tmp*2.0;
       value = (value-G*0.5)/1000.0;
       tmp = value - floor(value);
       B = tmp*2.0;
       value = (value-B*0.5)/1000.0;
       tmp = value - floor(value);
       A = tmp*2.0;
       color = new G4VisAttributes(G4Color(R,G,B,A));
       if (A==0.0) color->SetVisibility(false);
       (*lvciter)->SetVisAttributes(color);
     }
   }
   //
   // End of Auxiliary Information block
   //



 // return the world physical volume ----------------------------------------

 G4cout << G4endl << "The geometrical trees defined are : " << G4endl << G4endl;
 DumpGeometricalTree(fWorld,G4cout);

 return fWorld;
}


void SEMDetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aVolume,ostream& out,G4int depth)
{
  for(int isp=0;isp<depth;isp++)
  { out << "  "; }
  out << aVolume->GetName() << "[" << aVolume->GetCopyNo() << "] "
         << aVolume->GetLogicalVolume()->GetName() << " "
         << aVolume->GetLogicalVolume()->GetNoDaughters() << " "
         << aVolume->GetLogicalVolume()->GetMaterial()->GetName() << " "
         << "Atom density = " << aVolume->GetLogicalVolume()->GetMaterial()->GetTotNbOfAtomsPerVolume()*(cm*cm*cm) << " 1/cm3 ";
  if(aVolume->GetLogicalVolume()->GetSensitiveDetector())
  {
    out << " " << aVolume->GetLogicalVolume()->GetSensitiveDetector()
                            ->GetFullPathName();
  }
  G4Material* material = aVolume->GetLogicalVolume()->GetMaterial();
  G4MaterialPropertiesTable* aMPT = material->GetMaterialPropertiesTable();
  if(aMPT && material->GetState()==kStateGas) {
     out << " gas fractions";
     char property[10]="";
     G4int gas = 0;
     sprintf (property,"FRACTION%d",gas);
     while (aMPT->ConstPropertyExists(property)) {
        out << " " << aMPT->GetConstProperty(property);
        gas++;
        sprintf (property,"FRACTION%d",gas);
     }
  }
  out << G4endl;
  for(int i=0;i<aVolume->GetLogicalVolume()->GetNoDaughters();i++)
  { DumpGeometricalTree(aVolume->GetLogicalVolume()->GetDaughter(i),out,depth+1); }  
}


void SEMDetectorConstruction::DumpGDMLGeometry()
{
  G4VPhysicalVolume* g4wv =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking()->GetWorldVolume();

  parser->Write("SEM4.9.2_out.gdml",g4wv,"../GDMLSchema/gdml.xsd");

  //Old (4.9.0.p01):
  //G4GDMLWriter g4writer("../GDMLSchema/gdml.xsd", "SEMout.gdml",2);
  //g4writer.SetDetectorList(detectorlist);
  //g4writer.DumpGeometryInfo(g4wv);
}

G4LogicalVolume* SEMDetectorConstruction::FindLogicalVolume( const G4String& vn ) {

   return G4LogicalVolumeStore::GetInstance()->GetVolume(vn);

  // Old v4.9.0.p01:
  // //EK, GetLogicalVolume -> FindLogicalVolume (new method to deal with GDML modularity)
  //const G4LogicalVolume* lv = GDMLProcessor::GetInstance()->FindLogicalVolume( vn );
  //if (lv==0) G4cout << "Found zero logical volume " << vn << G4endl;
  //return const_cast<G4LogicalVolume*> (lv);
}

void SEMDetectorConstruction::setMaterial (G4String volumename, G4String materialname)
{
     // search the material by its name 
  G4Material* material = G4Material::GetMaterial(materialname);  
  if (material)
     {
       G4LogicalVolume* volume = FindLogicalVolume( volumename );
       if (volume) {
          volume->SetMaterial(material);
          G4cout << "Material of volume " << volumename << " set to " << materialname << G4endl;
          G4cout << "New geometrical tree: " << G4endl;
          DumpGeometricalTree(fWorld,G4cout);
       } else G4cerr << "SEMDetectorConstruction: failed to change material, volume " << volumename << " not found." << G4endl;
    }  else G4cerr << "SEMDetectorConstruction: failed to change material, material " << materialname << " not found." << G4endl;
}

void SEMDetectorConstruction::setFraction (G4String materialname, G4int index, G4double fraction) {
   G4Material* material = G4Material::GetMaterial(materialname);  
   if (material)
   {
      G4MaterialPropertiesTable* aMPT = material->GetMaterialPropertiesTable();
      if (!aMPT) {
         aMPT = new G4MaterialPropertiesTable();
         material->SetMaterialPropertiesTable(aMPT);
      }
      char property[10]="";
      sprintf (property,"FRACTION%d",index);
      aMPT->AddConstProperty(property,fraction);
      G4RunManagerKernel::GetRunManagerKernel()->PhysicsHasBeenModified();
   } else G4cerr << "SEMDetectorConstruction: Failed to change fraction, material " << materialname << " not found!" << G4endl;
}

void SEMDetectorConstruction::setGasDensity (G4String materialname, G4double density) {
   G4Material* material = G4Material::GetMaterial(materialname);  
   if (material)
   {
      G4MaterialPropertiesTable* aMPT = material->GetMaterialPropertiesTable();
      if (!aMPT) {
         aMPT = new G4MaterialPropertiesTable();
         material->SetMaterialPropertiesTable(aMPT);
      }
      aMPT->AddConstProperty("GASDENSITY",density);
      G4RunManagerKernel::GetRunManagerKernel()->PhysicsHasBeenModified();
   } else G4cerr << "SEMDetectorConstruction: Failed to change density, material " << materialname << " not found!" << G4endl;
}

G4double SEMDetectorConstruction::getGasDensity (G4String materialname) {
   G4Material* material = G4Material::GetMaterial(materialname);  
   if (material)
   {
      G4MaterialPropertiesTable* aMPT = material->GetMaterialPropertiesTable();
      if (aMPT) {
         if(aMPT->ConstPropertyExists("GASDENSITY")) {
            G4double density = aMPT->GetConstProperty("GASDENSITY");
            return density;
         }
      }
      return 0.;
   }
   return -1.;
}

void SEMDetectorConstruction::PrintTree (ostream& out)
{
   DumpGeometricalTree(fWorld, out);
}
