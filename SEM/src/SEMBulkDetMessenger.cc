// SEMBulkDetMessenger.cc
//

#include "SEMBulkDetMessenger.hh"
#include "SEMBulkDet.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

SEMBulkDetMessenger::SEMBulkDetMessenger(SEMBulkDet* Det)
:target(Det)
{
  G4String entry;
  entry="/detectors/"+Det->GetName()+"/";
  detDirectory = new G4UIdirectory(entry);
  detDirectory->SetGuidance("Commands specific for Det detector");

  entry="/detectors/"+Det->GetName()+"/multipleHits";
  detMultipleHitsCmd = new G4UIcmdWithABool(entry,this);
  detMultipleHitsCmd->SetGuidance("If true every hit will be reported.");
  detMultipleHitsCmd->SetGuidance("If false, only the first or last hit will be reported, depending on detectAtEnd,");
  detMultipleHitsCmd->SetGuidance("the deposited energy is summed over the whole track within the material.");
  detMultipleHitsCmd->SetParameterName("multiplehits",true);
  detMultipleHitsCmd->SetDefaultValue(false);

  entry="/detectors/"+Det->GetName()+"/detectAtEnd";
  detDetectAtEndCmd = new G4UIcmdWithABool(entry,this);
  detDetectAtEndCmd->SetGuidance("Turning this option on tells the detector to report the data at the end of each step.");
  detDetectAtEndCmd->SetParameterName("detectAtEnd",true);
  detDetectAtEndCmd->SetDefaultValue(false);

  entry="/detectors/"+Det->GetName()+"/detectTotal";
  detDetectTotalCmd = new G4UIcmdWithABool(entry,this);
  detDetectTotalCmd->SetGuidance("If true: Report only the total energy deposit.");
  detDetectTotalCmd->SetGuidance("It is reported at the position of the first hit in each event.");
  detDetectTotalCmd->SetGuidance("This overrules the multipleHits and detectAtEnd options.");
  detDetectTotalCmd->SetParameterName("detectTotal",true);
  detDetectTotalCmd->SetDefaultValue(false);

  entry="/detectors/"+Det->GetName()+"/energyCut";
  detEnergyCutCmd = new G4UIcmdWithADoubleAndUnit(entry,this);
  detEnergyCutCmd->SetGuidance("Set lower limit of deposited energy for a step to be stored when multipleHits=true.");
  detEnergyCutCmd->SetParameterName("Emin",false,false);
  detEnergyCutCmd->SetDefaultValue(0.0);
  detEnergyCutCmd->SetDefaultUnit("MeV");

  entry="/detectors/"+Det->GetName()+"/outputCounts";
  detCountsOutputCmd = new G4UIcmdWithAString(entry,this);
  detCountsOutputCmd->SetGuidance("Output number of hits on detector to given filename.");
  detCountsOutputCmd->SetGuidance("Enter a filename as parameter");
  detCountsOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Det->GetName()+"/outputHits";
  detHitsOutputCmd = new G4UIcmdWithAString(entry,this);
  detHitsOutputCmd->SetGuidance("Output hits of detector to given filename.");
  detHitsOutputCmd->SetGuidance("Enter a filename as parameter");
  detHitsOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Det->GetName()+"/outputEnergyHistogram";
  detHistOutputCmd = new G4UIcommand(entry,this);
  detHistOutputCmd->SetGuidance("Output energy histogram.");
  detHistOutputCmd->SetGuidance("First parameter: filename");
  detHistOutputCmd->SetGuidance("Second parameter: bin size in eV");
  detHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("filename",'s',false);
  detHistOutputCmd->SetParameter(p1);
  G4UIparameter* p2 = new G4UIparameter("binsize",'d',true);
  p2->SetDefaultValue("0.5");  // Default is 0.5 eV bins
  p2->SetParameterRange("binsize >= 0.0");
  detHistOutputCmd->SetParameter(p2);

  entry="/detectors/"+Det->GetName()+"/outputEnergyPDF";
  detEnergyPDFOutputCmd = new G4UIcmdWithAString(entry,this);
  detEnergyPDFOutputCmd->SetGuidance("Output probability density function for the energy of particles hitting the detector.");
  detEnergyPDFOutputCmd->SetGuidance("Enter a filename as parameter");
  detEnergyPDFOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Det->GetName()+"/outputRhoPDF";
  detRhoPDFOutputCmd = new G4UIcmdWithAString(entry,this);
  detRhoPDFOutputCmd->SetGuidance("Output of impact radius PDF.");
  detRhoPDFOutputCmd->SetGuidance("Enter a filename as parameter");
  detRhoPDFOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Det->GetName()+"/outputAngleHistogram";
  detAngleHistOutputCmd = new G4UIcommand(entry,this);
  detAngleHistOutputCmd->SetGuidance("Output histogram of momentum direction.");
  detAngleHistOutputCmd->SetGuidance("First parameter: filename");
  detAngleHistOutputCmd->SetGuidance("Second parameter (vector: 3 numbers): reference axis (no need to normalize)");
  detAngleHistOutputCmd->SetGuidance("Third parameter: number of bins (optional)");
  detAngleHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* q1 = new G4UIparameter("filename",'s',false);
  detAngleHistOutputCmd->SetParameter(q1);
  G4UIparameter* q2 = new G4UIparameter("x-component",'d',false);
  detAngleHistOutputCmd->SetParameter(q2);
  G4UIparameter* q3 = new G4UIparameter("y-component",'d',false);
  detAngleHistOutputCmd->SetParameter(q3);
  G4UIparameter* q4 = new G4UIparameter("z-component",'d',false);
  detAngleHistOutputCmd->SetParameter(q4);
  G4UIparameter* q5 = new G4UIparameter("number of bins",'i',true);
  q5->SetDefaultValue("20");  // Default is 20 bins
  detAngleHistOutputCmd->SetParameter(q5);

  entry="/detectors/"+Det->GetName()+"/outputDeposits";
  detDepositsOutputCmd = new G4UIcommand(entry,this);
  detDepositsOutputCmd->SetGuidance("Output energy deposits in detector");
  detDepositsOutputCmd->SetGuidance("First parameter: filename");
  G4UIparameter* r1 = new G4UIparameter("filename",'s',false);
  detDepositsOutputCmd->SetParameter(r1);

  entry="/detectors/"+Det->GetName()+"/outputTOFHistogram";
  detTimeHistOutputCmd = new G4UIcommand(entry,this);
  detTimeHistOutputCmd->SetGuidance("Output Time Of Flight histogram.");
  detTimeHistOutputCmd->SetGuidance("First parameter: filename");
  detTimeHistOutputCmd->SetGuidance("Second parameter: bin size in ns");
  detTimeHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p3 = new G4UIparameter("filename",'s',false);
  detTimeHistOutputCmd->SetParameter(p3);
  G4UIparameter* p4 = new G4UIparameter("binsize",'d',true);
  p4->SetDefaultValue("0.5");  // Default is 0.5 eV bins
  p4->SetParameterRange("binsize >= 0.0");
  detTimeHistOutputCmd->SetParameter(p4);

}

SEMBulkDetMessenger::~SEMBulkDetMessenger()
{
  delete detMultipleHitsCmd;
  delete detDetectAtEndCmd;
  delete detDetectTotalCmd;
  delete detEnergyCutCmd;
  delete detCountsOutputCmd;
  delete detHitsOutputCmd;
  delete detHistOutputCmd;
  delete detEnergyPDFOutputCmd;
  delete detRhoPDFOutputCmd;
  delete detAngleHistOutputCmd;
  delete detDepositsOutputCmd;
  delete detTimeHistOutputCmd;

  delete detDirectory;
}

void SEMBulkDetMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==detMultipleHitsCmd ) { target->SetMultipleHits(detMultipleHitsCmd->GetNewBoolValue(newValue)); }
  if( command==detDetectAtEndCmd ) { target->SetDetectAtEnd(detDetectAtEndCmd->GetNewBoolValue(newValue)); }
  if( command==detDetectTotalCmd ) { target->SetDetectTotal(detDetectTotalCmd->GetNewBoolValue(newValue)); }
  if( command==detEnergyCutCmd ) { target->SetEnergyCut(detEnergyCutCmd->GetNewDoubleValue(newValue)); }

  if( command==detCountsOutputCmd ) { target->OutputCounts(newValue); }
  if( command==detHitsOutputCmd ) { target->OutputHits(newValue); }
  if( command==detDepositsOutputCmd ) {
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename;
      target->OutputDeposits(filename);
  }
  if( command==detHistOutputCmd ) { 
      G4double binsize; 
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> binsize;
      target->OutputEnergyHistogram(filename,binsize); 
  }
  if( command==detEnergyPDFOutputCmd ) { target->OutputEnergyPDF(newValue); }
  if( command==detRhoPDFOutputCmd ) { target->OutputRhoPDF(newValue); }
  if( command==detAngleHistOutputCmd ) {
      string   filename;
      double   x,y,z;
      G4int    nbin;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> x >> y >> z >> nbin;
      G4ThreeVector normal(x,y,z);
      target->OutputAngleHistogram(filename,normal,nbin);
  }
  if( command==detTimeHistOutputCmd ) { 
      G4double binsize; 
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> binsize;
      target->OutputTimeHistogram(filename,binsize); 
  }

}

G4String SEMBulkDetMessenger::GetCurrentValue(G4UIcommand * /*command*/)
{
 G4String cv = "";
 return cv;
}
