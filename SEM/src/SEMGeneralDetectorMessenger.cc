// Jun 14 2007: Added commands that control the angular window selection

// SEMGeneralDetectorMessenger.cc
//

#include "SEMGeneralDetectorMessenger.hh"

#include "SEMGeneralDetector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"

SEMGeneralDetectorMessenger::SEMGeneralDetectorMessenger(SEMGeneralDetector * GD)
:target(GD)
{
  G4String entry;
  entry="/detectors/"+GD->GetName()+"/model";
  ModelCmd = new G4UIcmdWithAnInteger(entry,this);
  ModelCmd->SetGuidance("Select the model for conversion from hits to detector output (e.g. current)");
  ModelCmd->SetGuidance("Available models :");
  ModelCmd->SetGuidance("0: No model available");
  ModelCmd->SetGuidance("1: Solid state STEM-II detector");
  ModelCmd->SetParameterName("Model",true);
  ModelCmd->SetRange("Model>=0");
  ModelCmd->SetDefaultValue(0); // Default is that the model is not defined (i.e. no currents available)

  entry="/detectors/"+GD->GetName()+"/energyWindow";
  EwindowCmd = new G4UIcmdWithABool(entry,this);
  EwindowCmd->SetGuidance("Set an energy window for all relevant output commands for all detectors");
  EwindowCmd->SetParameterName("Ewindow",true);
  EwindowCmd->SetDefaultValue(false);

  entry="/detectors/"+GD->GetName()+"/minEnergy";
  minECmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  minECmd->SetGuidance("Set lower limit of kinetic energy for window.");
  minECmd->SetGuidance("A value smaller than 0 means a delta-energy w.r.t. the maximum energy.");
  minECmd->SetParameterName("Emin",true);
  minECmd->SetDefaultUnit("eV");
  minECmd->SetDefaultValue(0.0);

  entry="/detectors/"+GD->GetName()+"/maxEnergy";
  maxECmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  maxECmd->SetGuidance("Set upper limit of kinetic energy for window.");
  maxECmd->SetParameterName("Emax",true);
  maxECmd->SetDefaultUnit("eV");
  maxECmd->SetDefaultValue(1.0e30);
  maxECmd->SetRange("Emax>=0.0");

  entry="/detectors/"+GD->GetName()+"/angleWindow";
  AwindowCmd = new G4UIcmdWithABool(entry,this);
  AwindowCmd->SetGuidance("Set an angle window for all relevant output commands for all detectors");
  AwindowCmd->SetParameterName("Awindow",true);
  AwindowCmd->SetDefaultValue(false);

  entry="/detectors/"+GD->GetName()+"/minAngle";
  minACmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  minACmd->SetGuidance("Set lower limit of angle (degrees).");
  minACmd->SetParameterName("Amin",true);
  minACmd->SetDefaultUnit("deg");
  minACmd->SetDefaultValue(-180.0);

  entry="/detectors/"+GD->GetName()+"/maxAngle";
  maxACmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  maxACmd->SetGuidance("Set upper limit of angle.");
  maxACmd->SetParameterName("Amax",true);
  maxACmd->SetDefaultUnit("deg");
  maxACmd->SetDefaultValue(180.0);

  entry="/detectors/"+GD->GetName()+"/refAxis";
  dirAVecCmd = new G4UIcmdWith3Vector(entry,this);
  dirAVecCmd->SetGuidance("Reference axis for angular window");
  dirAVecCmd->SetParameterName("NX","NY","NZ",true,true);

  entry="/detectors/"+GD->GetName()+"/selectOutput";
  selectOutputCmd = new G4UIcmdWithAString(entry,this);
  selectOutputCmd->SetGuidance("Select output columns for outputHits command");
  selectOutputCmd->SetGuidance("Choices:");
  selectOutputCmd->SetGuidance("EID = Event ID");
  selectOutputCmd->SetGuidance("TID = Track ID");
  selectOutputCmd->SetGuidance("PID = Parent ID");
  selectOutputCmd->SetGuidance("Etot = Total energy");
  selectOutputCmd->SetGuidance("Ekin = Kinetic energy");
  selectOutputCmd->SetGuidance("Edeposit = Deposited energy");
  selectOutputCmd->SetGuidance("Eloss = Energy loss");
  selectOutputCmd->SetGuidance("Time = Time of hit");
  selectOutputCmd->SetGuidance("Pos = Position of hit");
  selectOutputCmd->SetGuidance("VPos = Origin of particle hitting the detector");
  selectOutputCmd->SetGuidance("Dir = Direction of particle hitting the detector");
  selectOutputCmd->SetGuidance("Depth = Minimum z-coordinate encountered for this particle");
  selectOutputCmd->SetGuidance("Maxr = Maximum distance from z-axis encountered for this particle");
  selectOutputCmd->SetGuidance("ScatterAngle = Scattering angle w.r.t. z-axis (between 0 and pi/2)");
  selectOutputCmd->SetGuidance("LV = Logical volume where the particle was created");
  selectOutputCmd->SetGuidance("Process = Process that created the particle");
  selectOutputCmd->SetGuidance("Category = Category of the particle (SE/BSE/gamma)");
  selectOutputCmd->SetGuidance("Potential = Potential at hit position");
  selectOutputCmd->SetGuidance("Default for outputHits: EID TID PID Etot Ekin Time Pos VPos Dir Depth Maxr LV Process Category");
  entry="/detectors/"+GD->GetName()+"/deselectOutput";
  deselectOutputCmd = new G4UIcmdWithAString(entry,this);
  deselectOutputCmd->SetGuidance("Deselect output columns for outputHits command");
  deselectOutputCmd->SetGuidance("Choices:");
  deselectOutputCmd->SetGuidance("EID = Event ID");
  deselectOutputCmd->SetGuidance("TID = Track ID");
  deselectOutputCmd->SetGuidance("PID = Parent ID");
  deselectOutputCmd->SetGuidance("Etot = Total energy");
  deselectOutputCmd->SetGuidance("Ekin = Kinetic energy");
  deselectOutputCmd->SetGuidance("Edeposit = Deposited energy");
  deselectOutputCmd->SetGuidance("Eloss = Energy loss");
  deselectOutputCmd->SetGuidance("Time = Time of hit");
  deselectOutputCmd->SetGuidance("Pos = Position of hit");
  deselectOutputCmd->SetGuidance("VPos = Origin of particle hitting the detector");
  deselectOutputCmd->SetGuidance("Dir = Direction of particle hitting the detector");
  deselectOutputCmd->SetGuidance("Depth = Minimum z-coordinate encountered for this particle");
  deselectOutputCmd->SetGuidance("Maxr = Maximum distance from z-axis encountered for this particle");
  deselectOutputCmd->SetGuidance("ScatterAngle = Scattering angle w.r.t. z-axis (between 0 and pi/2)");
  deselectOutputCmd->SetGuidance("LV = Logical volume where the particle was created");
  deselectOutputCmd->SetGuidance("Process = Process that created the particle");
  deselectOutputCmd->SetGuidance("Category = Category of the particle (SE/BSE/gamma)");
  deselectOutputCmd->SetGuidance("Potential = Potential at hit position");
  deselectOutputCmd->SetGuidance("All = Unselects all columns");
}

SEMGeneralDetectorMessenger::~SEMGeneralDetectorMessenger()
{
  delete ModelCmd;
  delete EwindowCmd;
  delete minECmd;
  delete maxECmd;
  delete AwindowCmd;
  delete minACmd;
  delete maxACmd;
  delete dirAVecCmd;
  delete selectOutputCmd;
  delete deselectOutputCmd;
}

void SEMGeneralDetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==ModelCmd ) { target->SetDetectorModel(ModelCmd->GetNewIntValue(newValue)); }
  if( command==EwindowCmd ) { target->SetEwindow(EwindowCmd->GetNewBoolValue(newValue)); }
  if( command==minECmd ) { target->SetEmin(minECmd->GetNewDoubleValue(newValue)); }
  if( command==maxECmd ) { target->SetEmax(maxECmd->GetNewDoubleValue(newValue)); }
  if( command==AwindowCmd ) { target->SetAwindow(AwindowCmd->GetNewBoolValue(newValue)); }
  if( command==minACmd ) { target->SetAmin(minACmd->GetNewDoubleValue(newValue)); }
  if( command==maxACmd ) { target->SetAmax(maxACmd->GetNewDoubleValue(newValue)); }
  if( command==dirAVecCmd ) { target->SetAdir(dirAVecCmd->GetNew3VectorValue(newValue)); }
  if( command==selectOutputCmd ) { 
     G4String materialname;
     const char* nv = (const char*)newValue;
     std::istringstream is(nv);
     G4String Selection;
     while(!is.eof()) {
        is >> Selection;
        if (Selection=="EID") target->SetOutputEID();
        if (Selection=="TID") target->SetOutputTID();
        if (Selection=="PID") target->SetOutputPID();
        if (Selection=="Etot") target->SetOutputEtot();
        if (Selection=="Ekin") target->SetOutputEkin();
        if (Selection=="Edeposit") target->SetOutputEdeposit();
        if (Selection=="Eloss") target->SetOutputEloss();
        if (Selection=="Time") target->SetOutputTime();
        if (Selection=="Pos") target->SetOutputPosition();
        if (Selection=="VPos") target->SetOutputVertexPosition();
        if (Selection=="Dir") target->SetOutputDirection();
        if (Selection=="Depth") target->SetOutputDepth();
        if (Selection=="Maxr") target->SetOutputMaxRadius();
        if (Selection=="r") target->SetOutputRadius();
        if (Selection=="Potential") target->SetOutputPotential();
        if (Selection=="ScatterAngle") target->SetOutputScatteringAngle();
        if (Selection=="LV") target->SetOutputLogicalVolume();
        if (Selection=="Process") target->SetOutputProcess();
        if (Selection=="Category") target->SetOutputCategory();
     }
  }
  if( command==deselectOutputCmd ) { 
     G4String materialname;
     const char* nv = (const char*)newValue;
     std::istringstream is(nv);
     G4String Selection;
     while(!is.eof()) {
        is >> Selection;
        if (Selection=="EID") target->UnsetOutputEID();
        if (Selection=="TID") target->UnsetOutputTID();
        if (Selection=="PID") target->UnsetOutputPID();
        if (Selection=="Etot") target->UnsetOutputEtot();
        if (Selection=="Ekin") target->UnsetOutputEkin();
        if (Selection=="Edeposit") target->UnsetOutputEdeposit();
        if (Selection=="Eloss") target->UnsetOutputEloss();
        if (Selection=="Time") target->UnsetOutputTime();
        if (Selection=="Pos") target->UnsetOutputPosition();
        if (Selection=="VPos") target->UnsetOutputVertexPosition();
        if (Selection=="Dir") target->UnsetOutputDirection();
        if (Selection=="Depth") target->UnsetOutputDepth();
        if (Selection=="Maxr") target->UnsetOutputMaxRadius();
        if (Selection=="r") target->UnsetOutputRadius();
        if (Selection=="Potential") target->UnsetOutputPotential();
        if (Selection=="ScatterAngle") target->UnsetOutputScatteringAngle();
        if (Selection=="LV") target->UnsetOutputLogicalVolume();
        if (Selection=="Process") target->UnsetOutputProcess();
        if (Selection=="Category") target->UnsetOutputCategory();
        if (Selection=="All") {
           target->UnsetOutputEID();
           target->UnsetOutputTID();
           target->UnsetOutputPID();
           target->UnsetOutputEtot();
           target->UnsetOutputEkin();
           target->UnsetOutputEdeposit();
           target->UnsetOutputEloss();
           target->UnsetOutputTime();
           target->UnsetOutputPosition();
           target->UnsetOutputVertexPosition();
           target->UnsetOutputDirection();
           target->UnsetOutputDepth();
           target->UnsetOutputMaxRadius();
           target->UnsetOutputRadius();
           target->UnsetOutputPotential();
           target->UnsetOutputScatteringAngle();
           target->UnsetOutputLogicalVolume();
           target->UnsetOutputProcess();
           target->UnsetOutputCategory();
        }
     }
  }
}

G4String SEMGeneralDetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==ModelCmd ) { cv = ModelCmd->ConvertToString(target->GetDetectorModel()); }
  if( command==EwindowCmd ) { cv = EwindowCmd->ConvertToString(target->GetEwindow()); }
  if( command==minECmd ) { cv = minECmd->ConvertToString(target->GetEmin()); }
  if( command==maxECmd ) { cv = maxECmd->ConvertToString(target->GetEmax()); }
  if( command==AwindowCmd ) { cv = AwindowCmd->ConvertToString(target->GetAwindow()); }
  if( command==minACmd ) { cv = minACmd->ConvertToString(target->GetAmin()); }
  if( command==maxACmd ) { cv = maxACmd->ConvertToString(target->GetAmax()); }

  return cv;
}
