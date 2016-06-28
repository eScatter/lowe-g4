// SEMCounterMessenger.cc
//

// Oct 16 2008 Created

#include "SEMCounterMessenger.hh"
#include "SEMCounter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"

SEMCounterMessenger::SEMCounterMessenger(SEMCounter* Counter)
:target(Counter)
{
  G4String entry;
  entry="/detectors/"+Counter->GetName()+"/";
  counterDirectory = new G4UIdirectory(entry);
  entry="Commands specific for "+Counter->GetName()+" detector";
  counterDirectory->SetGuidance(entry);

  entry="/detectors/"+Counter->GetName()+"/outputCounts";
  counterCountsOutputCmd = new G4UIcmdWithAString(entry,this);
  counterCountsOutputCmd->SetGuidance("Output number of hits on detector to given filename.");
  counterCountsOutputCmd->SetGuidance("Enter a filename as parameter");
  counterCountsOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Counter->GetName()+"/outputHits";
  counterHitsOutputCmd = new G4UIcmdWithAString(entry,this);
  counterHitsOutputCmd->SetGuidance("Output hits of detector to given filename.");
  counterHitsOutputCmd->SetGuidance("Enter a filename as parameter");
  counterHitsOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Counter->GetName()+"/outputEnergyHistogram";
  counterHistOutputCmd = new G4UIcommand(entry,this);
  counterHistOutputCmd->SetGuidance("Output energy histogram.");
  counterHistOutputCmd->SetGuidance("First parameter: filename");
  counterHistOutputCmd->SetGuidance("Second parameter: bin size in eV");
  counterHistOutputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("filename",'s',false);
  counterHistOutputCmd->SetParameter(p1);
  G4UIparameter* p2 = new G4UIparameter("binsize",'d',true);
  p2->SetDefaultValue("0.5");  // Default is 0.5 eV bins
  p2->SetParameterRange("binsize >= 0.0");
  counterHistOutputCmd->SetParameter(p2);

  entry="/detectors/"+Counter->GetName()+"/outputEnergyPDF";
  counterEnergyPDFOutputCmd = new G4UIcmdWithAString(entry,this);
  counterEnergyPDFOutputCmd->SetGuidance("Output probability density function for the energy of particles hitting the detector.");
  counterEnergyPDFOutputCmd->SetGuidance("Enter a filename as parameter");
  counterEnergyPDFOutputCmd->SetParameterName("filename",false);

  entry="/detectors/"+Counter->GetName()+"/verbose";
  counterVerboseCmd = new G4UIcmdWithAnInteger(entry,this);
  counterVerboseCmd->SetGuidance("Verbose level for the counter");
  counterVerboseCmd->SetParameterName("level",true);
  counterVerboseCmd->SetRange("level>=0");
  counterVerboseCmd->SetDefaultValue(0);

}

SEMCounterMessenger::~SEMCounterMessenger()
{
  delete counterVerboseCmd;
  delete counterCountsOutputCmd;
  delete counterHitsOutputCmd;
  delete counterHistOutputCmd;
  delete counterEnergyPDFOutputCmd;
  delete counterDirectory;
}

void SEMCounterMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==counterVerboseCmd ) { target->SetVerbose(counterVerboseCmd->GetNewIntValue(newValue)); }
  if( command==counterCountsOutputCmd ) { target->OutputCounts(newValue); }
  if( command==counterHitsOutputCmd ) { target->OutputHits(newValue); }
  if( command==counterHistOutputCmd ) { 
      G4double binsize; 
      string   filename;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> filename >> binsize;
      target->OutputEnergyHistogram(filename,binsize); 
  }
  if( command==counterEnergyPDFOutputCmd ) { target->OutputEnergyPDF(newValue); }
}

G4String SEMCounterMessenger::GetCurrentValue(G4UIcommand * /*command*/)
{
  G4String cv;

  return cv;
}
