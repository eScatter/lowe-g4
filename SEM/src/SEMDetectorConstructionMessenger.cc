#include "SEMDetectorConstructionMessenger.hh"

#include "SEMDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

SEMDetectorConstructionMessenger::SEMDetectorConstructionMessenger(SEMDetectorConstruction* mpga)
:target(mpga)
{
  detectorsDirectory = new G4UIdirectory("/detectors/");
  detectorsDirectory->SetGuidance("Miscellaneous commands for SEM");

  writeCmd = new G4UIcommand("/detectors/writegdml",this);
  writeCmd->SetGuidance("Write current geometry to gdml file");

  materialCmd = new G4UIcommand("/detectors/setmaterial",this);
  materialCmd->SetGuidance("Select material of a volume.");
  materialCmd->SetGuidance("First parameter: name of the volume.");
  materialCmd->SetGuidance("Second parameter: name of the new material.");
  materialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("volumename",'s',false);
  materialCmd->SetParameter(p1);
  G4UIparameter* p2 = new G4UIparameter("materialname",'s',false);
  materialCmd->SetParameter(p2);

  gasfractionsCmd = new G4UIcmdWithAString("/detectors/gasfractions",this);
  gasfractionsCmd->SetGuidance("Set gas composition fractions for a material.");
  gasfractionsCmd->SetGuidance("First parameter: material name.");
  gasfractionsCmd->SetGuidance("Following parameters: relative gas fractions of the various gases.");

  gasdensityCmd = new G4UIcommand("/detectors/setgasdensity",this);
  gasdensityCmd->SetGuidance("Set number density (in units of mm^-3) of a gas.");
  gasdensityCmd->SetGuidance("First parameter: material name.");
  gasdensityCmd->SetGuidance("Second parameter: density (without unit).");
  G4UIparameter* p3 = new G4UIparameter("materialname",'s',false);
  gasdensityCmd->SetParameter(p3);
  G4UIparameter* p4 = new G4UIparameter("density",'d',true);
  gasdensityCmd->SetParameter(p4);
  p4->SetDefaultValue("0.0");
  p4->SetParameterRange("density >= 0.0");

  getdensityCmd = new G4UIcmdWithAString("/detectors/getgasdensity",this);
  getdensityCmd->SetGuidance("Get number density for a gas (in units of mm^-3).");
  getdensityCmd->SetGuidance("Parameter: material name.");
  
  printtreeCmd = new G4UIcommand("/detectors/printtree",this);
  printtreeCmd->SetGuidance("Display current geometrical tree(s)");
}

SEMDetectorConstructionMessenger::~SEMDetectorConstructionMessenger()
{
  delete printtreeCmd;
  delete materialCmd;
  delete gasfractionsCmd;
  delete gasdensityCmd;
  delete getdensityCmd;
  delete writeCmd;
  delete detectorsDirectory;
}

void SEMDetectorConstructionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
 if( command==writeCmd )
 { target->DumpGDMLGeometry(); } else
 if( command==materialCmd )
 {
	 G4String volumename; 
	 G4String materialname;
	 const char* nv = (const char*)newValue;
	 std::istringstream is(nv);
	 is >> volumename >> materialname;
	 //G4String materialname = "Gold";
	 target->setMaterial(volumename,materialname);
 } else
 if( command==gasfractionsCmd )
 {
	 G4String materialname;
	 const char* nv = (const char*)newValue;
	 std::istringstream is(nv);
	 is >> materialname;
	 G4double fraction;
	 G4int index = 0;
	 while(!is.eof()) {
		 is >> fraction;
		 target->setFraction(materialname,index,fraction);
		 index++;
	 }
	 G4cout << "New geometrical tree: " << G4endl;
	 target->PrintTree(G4cout);
 } else
 if( command==gasdensityCmd )
 {
	 G4String materialname;
	 G4double density;
	 const char* nv = (const char*)newValue;
	 std::istringstream is(nv);
	 is >> materialname >> density;
	 target->setGasDensity(materialname,density);
 } else
 if( command==getdensityCmd ) {
	 G4String materialname;
	 const char* nv = (const char*)newValue;
	 std::istringstream is(nv);
	 is >> materialname ;
	 G4double density = target->getGasDensity(materialname);
	 if (density >= 0.) {
		 G4cout << "Molecular number density of gas " << materialname << " is " << target->getGasDensity(materialname) << G4endl;
	 } else G4cout << "Error: material not found." << G4endl;
  } else
 if( command==printtreeCmd )
 {
	 target->PrintTree(G4cout);
 }
}

G4String SEMDetectorConstructionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;

  if( command==writeCmd ) { G4cout << "Not implemented!" << G4endl; }

  return cv;
}

