// SEMPlaneMessenger.cc
//

// Feb-28 2012 Added second way of input for direction of detector: theta and phi angles
// Sep-27 2010 Added output of energy deposits
// May-21 2007 Added energy selection for Angle histogram, default is no energy selection

#include "SEMPlaneMessenger.hh"
#include "SEMPlane.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

SEMPlaneMessenger::SEMPlaneMessenger(SEMPlane* Plane)
:target(Plane)
{
  G4String entry;
  entry="/detectors/"+Plane->GetName()+"/";
  planeDirectory = new G4UIdirectory(entry);
  entry="Commands specific for "+Plane->GetName()+" detector";
  planeDirectory->SetGuidance(entry);

  entry="/detectors/"+Plane->GetName()+"/multipleHits";
  multiCmd = new G4UIcmdWithABool(entry,this);
  multiCmd->SetGuidance("Determines if detector will record multiple hits of the same particle.");
  multiCmd->SetParameterName("multiplehits",true);
  multiCmd->SetDefaultValue(false);

  entry="/detectors/"+Plane->GetName()+"/detectStopped";
  pairCmd = new G4UIcmdWithABool(entry,this);
  pairCmd->SetGuidance("Turning this option on tells the detector to detect only particles that end their track in the detector.");
  pairCmd->SetGuidance("Particles leaving the detector are not counted.");
  pairCmd->SetGuidance("The hit location is the position where the particle has lost all its energy and is no longer able");
  pairCmd->SetGuidance("to generate new secondaries.");
  pairCmd->SetParameterName("detectStopped",true);
  pairCmd->SetDefaultValue(false);

  entry="/detectors/"+Plane->GetName()+"/directional";
  dirCmd = new G4UIcmdWithABool(entry,this);
  dirCmd->SetGuidance("Determines if detector will be sensitive in a given direction.");
  dirCmd->SetParameterName("directional",true);
  dirCmd->SetDefaultValue(true);

  entry="/detectors/"+Plane->GetName()+"/direction";
  dirVecCmd = new G4UIcmdWith3Vector(entry,this);
  dirVecCmd->SetGuidance("Direction in which the detector is sensitive.");
  dirVecCmd->SetParameterName("X","Y","Z",true,true);

  dirTheta=0.0;
  entry="/detectors/"+Plane->GetName()+"/directionTheta";
  dirThetaCmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  dirThetaCmd->SetGuidance("Inclination angle theta for the direction in which the detector is sensitive");
  dirThetaCmd->SetParameterName("dirTheta",true,true);
  dirThetaCmd->SetDefaultUnit("deg");
  dirThetaCmd->SetUnitCandidates("rad deg");
  dirThetaCmd->SetDefaultValue(0.0);

  dirPhi=0.0;
  entry="/detectors/"+Plane->GetName()+"/directionPhi";
  dirPhiCmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  dirPhiCmd->SetGuidance("Azimuthal angle phi for the direction in which the detector is sensitive");
  dirPhiCmd->SetParameterName("dirPhi",true,true);
  dirPhiCmd->SetDefaultUnit("deg");
  dirPhiCmd->SetUnitCandidates("rad deg");
  dirPhiCmd->SetDefaultValue(0.0);

  entry="/detectors/"+Plane->GetName()+"/outputCounts";
  planeCountsOutputCmd = new G4UIcmdWithAString(entry,this);
  planeCountsOutputCmd->SetGuidance("Output number of hits on detector to given filename.");
  planeCountsOutputCmd->SetGuidance("Enter a filename as parameter");
  planeCountsOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Plane->GetName()+"/outputHits";
  planeHitsOutputCmd = new G4UIcmdWithAString(entry,this);
  planeHitsOutputCmd->SetGuidance("Output hits of detector to given filename.");
  planeHitsOutputCmd->SetGuidance("Enter a filename as parameter");
  planeHitsOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Plane->GetName()+"/outputEnergyHistogram";
  planeHistOutputCmd = new G4UIcommand(entry,this);
  planeHistOutputCmd->SetGuidance("Output energy histogram.");
  planeHistOutputCmd->SetGuidance("First parameter: filename");
  planeHistOutputCmd->SetGuidance("Second parameter: bin size in eV");
  planeHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("filename",'s',false);
  planeHistOutputCmd->SetParameter(p1);
  G4UIparameter* p2 = new G4UIparameter("binsize",'d',true);
  p2->SetDefaultValue("0.5");  // Default is 0.5 eV bins
  p2->SetParameterRange("binsize >= 0.0");
  planeHistOutputCmd->SetParameter(p2);

  entry="/detectors/"+Plane->GetName()+"/outputEnergyPDF";
  planeEnergyPDFOutputCmd = new G4UIcmdWithAString(entry,this);
  planeEnergyPDFOutputCmd->SetGuidance("Output probability density function for the energy of particles hitting the detector.");
  planeEnergyPDFOutputCmd->SetGuidance("Enter a filename as parameter");
  planeEnergyPDFOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Plane->GetName()+"/outputGammaPDF";
  planeGammaPDFOutputCmd = new G4UIcmdWithAString(entry,this);
  planeGammaPDFOutputCmd->SetGuidance("Output probability density function for the energy of photons hitting the detector.");
  planeGammaPDFOutputCmd->SetGuidance("Enter a filename as parameter");
  planeGammaPDFOutputCmd->SetParameterName("filename",false);


  entry="/detectors/"+Plane->GetName()+"/outputRhoHistogram";
  planeRhoHistOutputCmd = new G4UIcommand(entry,this);
  planeRhoHistOutputCmd->SetGuidance("Output of impact radius histogram.");
  planeRhoHistOutputCmd->SetGuidance("First parameter: filename");
  planeRhoHistOutputCmd->SetGuidance("Second parameter: bin size in nanometer");
  planeRhoHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* a1 = new G4UIparameter("filename",'s',false);
  planeRhoHistOutputCmd->SetParameter(a1);
  G4UIparameter* a2 = new G4UIparameter("binsize",'d',true);
  a2->SetDefaultValue("1.0");  // Default is 1 nm bins
  a2->SetParameterRange("binsize >= 0.0");
  planeRhoHistOutputCmd->SetParameter(a2);

  entry="/detectors/"+Plane->GetName()+"/outputRhoPDF";
  planeRhoPDFOutputCmd = new G4UIcmdWithAString(entry,this);
  planeRhoPDFOutputCmd->SetGuidance("Output of impact radius PDF.");
  planeRhoPDFOutputCmd->SetGuidance("Enter a filename as parameter");
  planeRhoPDFOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Plane->GetName()+"/outputAngleHistogram";
  planeAngleHistOutputCmd = new G4UIcommand(entry,this);
  planeAngleHistOutputCmd->SetGuidance("Output histogram of momentum direction.");
  planeAngleHistOutputCmd->SetGuidance("First parameter: filename");
  planeAngleHistOutputCmd->SetGuidance("Second parameter (vector: 3 numbers): reference axis (no need to normalize)");
  planeAngleHistOutputCmd->SetGuidance("Third parameter: number of bins (optional)");
  planeAngleHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* q1 = new G4UIparameter("filename",'s',false);
  planeAngleHistOutputCmd->SetParameter(q1);
  G4UIparameter* q2 = new G4UIparameter("x-component",'d',false);
  planeAngleHistOutputCmd->SetParameter(q2);
  G4UIparameter* q3 = new G4UIparameter("y-component",'d',false);
  planeAngleHistOutputCmd->SetParameter(q3);
  G4UIparameter* q4 = new G4UIparameter("z-component",'d',false);
  planeAngleHistOutputCmd->SetParameter(q4);
  G4UIparameter* q5 = new G4UIparameter("number of bins",'i',true);
  q5->SetDefaultValue("20");  // Default is 20 bins
  planeAngleHistOutputCmd->SetParameter(q5);

  entry="/detectors/"+Plane->GetName()+"/outputDeposits";
  planeDepositsOutputCmd = new G4UIcommand(entry,this);
  planeDepositsOutputCmd->SetGuidance("Output energy deposits in detector to given filename.");
  planeDepositsOutputCmd->SetGuidance("First parameter: filename");
  planeDepositsOutputCmd->SetGuidance("Second parameter: minimum energy deposit to be reported (eV)");
  G4UIparameter* r1 = new G4UIparameter("filename",'s',false);
  planeDepositsOutputCmd->SetParameter(r1);
  G4UIparameter* r2 = new G4UIparameter("mindeposit",'d',false);
  r2->SetDefaultValue("0.0");
  r2->SetParameterRange("mindeposit >= 0.0");
  planeDepositsOutputCmd->SetParameter(r2);

}

SEMPlaneMessenger::~SEMPlaneMessenger()
{
  delete multiCmd;
  delete pairCmd;
  delete dirCmd;
  delete dirVecCmd;
  delete planeCountsOutputCmd;
  delete planeHitsOutputCmd;
  delete planeDepositsOutputCmd;
  delete planeHistOutputCmd;
  delete planeEnergyPDFOutputCmd;
  delete planeGammaPDFOutputCmd;
  delete planeRhoHistOutputCmd;
  delete planeRhoPDFOutputCmd;
  delete planeAngleHistOutputCmd;
  delete dirThetaCmd;
  delete dirPhiCmd;

  delete planeDirectory;
}

void SEMPlaneMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==multiCmd ) { target->SetMultipleHits(multiCmd->GetNewBoolValue(newValue)); }
  if( command==pairCmd ) { target->SetCountEHPairs(pairCmd->GetNewBoolValue(newValue)); }
  if( command==dirCmd ) { target->SetDirectional(dirCmd->GetNewBoolValue(newValue)); }
  if( command==dirVecCmd ) {target->SetDirection(dirVecCmd->GetNew3VectorValue(newValue)); }
  if( command==dirThetaCmd ) {dirTheta=dirThetaCmd->GetNewDoubleValue(newValue); target->SetDirection(dirTheta,dirPhi); }
  if( command==dirPhiCmd ) {dirPhi=dirPhiCmd->GetNewDoubleValue(newValue); target->SetDirection(dirTheta,dirPhi); }
  if( command==planeCountsOutputCmd ) { target->OutputCounts(newValue); }
  if( command==planeHitsOutputCmd ) { target->OutputHits(newValue); }
  if( command==planeDepositsOutputCmd ) { 
      G4double mindeposit; 
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> mindeposit;
      mindeposit *= eV;
      target->OutputDeposits(filename,mindeposit); 
  }
  if( command==planeHistOutputCmd ) { 
      G4double binsize; 
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> binsize;
      target->OutputEnergyHistogram(filename,binsize); 
  }
  if( command==planeRhoHistOutputCmd ) { 
      G4double binsize; 
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> binsize;
      target->OutputRhoHistogram(filename,binsize); 
  }
  if( command==planeEnergyPDFOutputCmd ) { target->OutputEnergyPDF(newValue); }
  if( command==planeGammaPDFOutputCmd ) { target->OutputGammaPDF(newValue); }
  if( command==planeRhoPDFOutputCmd ) { target->OutputRhoPDF(newValue); }
  if( command==planeAngleHistOutputCmd ) { 
      string   filename;
      double   x,y,z;
      G4int    nbin;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> x >> y >> z >> nbin;
      G4ThreeVector normal(x,y,z);
      target->OutputAngleHistogram(filename,normal,nbin); 
  }

}

G4String SEMPlaneMessenger::GetCurrentValue(G4UIcommand * /*command*/)
{
  G4String cv;

  return cv;
}
