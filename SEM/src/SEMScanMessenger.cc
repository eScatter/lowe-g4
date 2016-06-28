// SEMScanMessenger.cc
//
#include "SEMScanMessenger.hh"

#include "SEMRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

#include "G4UIparameter.hh"

SEMScanMessenger::SEMScanMessenger(SEMRunAction* runAction)
:scanManager(runAction)
{
  scanDirectory = new G4UIdirectory("/scan/");
  scanDirectory->SetGuidance("Scan control commands.");
  
  scanCmd = new G4UIcommand("/scan/scanOn",this);
  scanCmd->SetGuidance("Start a Scan along a line.");
  scanCmd->AvailableForStates(G4State_Idle);
  G4UIparameter* s1 = new G4UIparameter("numOfEventPerPixel",'i',true);
  s1->SetDefaultValue(1);
  s1->SetParameterRange("numOfEventPerPixel >= 0");
  scanCmd->SetParameter(s1);
  G4UIparameter* s2 = new G4UIparameter("filename",'s',true);
  s2->SetDefaultValue("");
  scanCmd->SetParameter(s2);

  imageCmd = new G4UIcommand("/scan/imageOn",this);
  imageCmd->SetGuidance("Start an Image scan.");
  imageCmd->SetGuidance("- n: number of primaries per pixel");
  imageCmd->SetGuidance("- filename: output filename for intermediate results");
  imageCmd->AvailableForStates(G4State_Idle);
  G4UIparameter* t1 = new G4UIparameter("numOfEventPerPixel",'i',true);
  t1->SetDefaultValue(1);
  t1->SetParameterRange("numOfEventPerPixel >= 0");
  imageCmd->SetParameter(t1);
  G4UIparameter* tt1 = new G4UIparameter("filename",'s',true);
  tt1->SetDefaultValue("");
  imageCmd->SetParameter(tt1);

  scanXdirCmd = new G4UIcmdWith3Vector("/scan/scanXdir",this);
  scanXdirCmd->SetGuidance("Total vector for 1st scan direction (X-direction of image).");
  scanXdirCmd->SetParameterName("X","Y","Z",true,true);

  scanYdirCmd = new G4UIcmdWith3Vector("/scan/scanYdir",this);
  scanYdirCmd->SetGuidance("Total vector for 2nd scan direction (Y-direction of image.");
  scanYdirCmd->SetParameterName("X","Y","Z",true,true);

  scanNxCmd = new G4UIcmdWithAnInteger("/scan/scanNx",this);
  scanNxCmd->SetGuidance("Number of points in the 1st scan direction (X-direction of image).");
  scanNxCmd->SetParameterName("Nx",true,true);
  scanNxCmd->SetRange("Nx>0");
  scanNxCmd->SetDefaultValue(1);
  
  scanNyCmd = new G4UIcmdWithAnInteger("/scan/scanNy",this);
  scanNyCmd->SetGuidance("Number of points in the 2nd scan direction (Y-direction of image).");
  scanNyCmd->SetParameterName("Ny",true,true);
  scanNyCmd->SetRange("Ny>0");
  scanNyCmd->SetDefaultValue(1);

  scanLxCmd = new G4UIcmdWithADoubleAndUnit("/scan/scanLx",this);
  scanLxCmd->SetGuidance("Length of image in the 1st scan direction (X-direction of image).");
  scanLxCmd->SetParameterName("Lx",true,true);
  scanLxCmd->SetRange("Lx>0");
  scanLxCmd->SetDefaultUnit("mm");
  scanLxCmd->SetDefaultValue(0.0);

  scanLyCmd = new G4UIcmdWithADoubleAndUnit("/scan/scanLy",this);
  scanLyCmd->SetGuidance("Length of image in the 2nd scan direction (Y-direction of image).");
  scanLyCmd->SetParameterName("Ly",true,true);
  scanLyCmd->SetRange("Ly>0");
  scanLyCmd->SetDefaultUnit("mm");
  scanLyCmd->SetDefaultValue(0.0);

  scanOriginCmd = new G4UIcmdWith3VectorAndUnit("/scan/scanOrigin",this);
  scanOriginCmd->SetGuidance("Starting position for scan or image.");
  scanOriginCmd->SetParameterName("X","Y","Z",true,true);
  scanOriginCmd->SetDefaultUnit("nm");

  scanOutputCmd = new G4UIcommand("/scan/outputLinescan",this);
  scanOutputCmd->SetGuidance("Output of linescan to a file");
  scanOutputCmd->SetGuidance("Takes a filename and the name of one of your detectors.");
  scanOutputCmd->AvailableForStates(G4State_Idle);
  G4UIparameter* t2 = new G4UIparameter("filename",'s',false);
  scanOutputCmd->SetParameter(t2);
  G4UIparameter* t3 = new G4UIparameter("detector name",'s',false);
  scanOutputCmd->SetParameter(t3);

  imageOutputCmd = new G4UIcommand("/scan/outputImage",this);
  imageOutputCmd->SetGuidance("Output of Image to a pgm file");
  imageOutputCmd->SetGuidance("Arguments:");
  imageOutputCmd->SetGuidance(" filename       :    filename for output");
  imageOutputCmd->SetGuidance(" detectorname   :    name of detector selected for output");
  imageOutputCmd->SetGuidance(" category       :    1 = Total number of electrons");
  imageOutputCmd->SetGuidance("                     2 = BSE");
  imageOutputCmd->SetGuidance("                     3 = SE");
  imageOutputCmd->SetGuidance("                     4 = SE3 (coupled detector only)");
  imageOutputCmd->SetGuidance("                     5 = Counts in Etot Window ");
  imageOutputCmd->SetGuidance("                     6 = Counts in Ekin Window ");
  imageOutputCmd->SetGuidance("                     7 = Counts in Angle Window ");
  imageOutputCmd->SetGuidance("                     8 = Current (only if a model is given)");
  imageOutputCmd->SetGuidance(" autoscale      :    true if you want autoscaling");
  imageOutputCmd->SetGuidance(" min            :    value that corresponds to black");
  imageOutputCmd->SetGuidance(" max            :    value that corresponds to white");
  imageOutputCmd->AvailableForStates(G4State_Idle);
  G4UIparameter* t4 = new G4UIparameter("filename",'s',false);
  imageOutputCmd->SetParameter(t4);

  G4UIparameter* t5 = new G4UIparameter("detector name",'s',false);
  imageOutputCmd->SetParameter(t5);

  G4UIparameter* t6 = new G4UIparameter("category",'i',true);
  t6->SetGuidance("1: nTotal  2: nBSE   3: nSE   4: nSE3  5: nEtotWindow  6: nEkinWindow  7: nAngWindow  8: Current");
  t6->SetDefaultValue(1);
  t6->SetParameterRange("(category>0) && (category<9)");
  imageOutputCmd->SetParameter(t6);

  G4UIparameter* t7 = new G4UIparameter("autoscale",'b',true);
  t7->SetGuidance("If true, the image is scaled automatically. If false, the user has to supply the minimum and maximum");
  t7->SetDefaultValue(true);
  imageOutputCmd->SetParameter(t7);

  G4UIparameter* t8 = new G4UIparameter("min",'i',true);
  t8->SetGuidance("Number of counts that will correspond to black in the pgm image");
  t8->SetDefaultValue(0.0);
  t8->SetParameterRange("min>=0.0");
  imageOutputCmd->SetParameter(t8);

  G4UIparameter* t9 = new G4UIparameter("max",'i',true);
  t9->SetGuidance("Number of counts that will correspond to white in the pgm image");
  t9->SetDefaultValue(0.0);
  t9->SetParameterRange("max>0.0");
  imageOutputCmd->SetParameter(t9);

  imageASCIIOutputCmd = new G4UIcommand("/scan/outputImageASCII",this);
  imageASCIIOutputCmd->SetGuidance("Output of Image to an ASCII file");
  imageASCIIOutputCmd->SetGuidance("Arguments:");
  imageASCIIOutputCmd->SetGuidance(" filename       :    filename for output");
  imageASCIIOutputCmd->SetGuidance(" detectorname   :    name of detector selected for output");
  imageASCIIOutputCmd->AvailableForStates(G4State_Idle);
  G4UIparameter* t10 = new G4UIparameter("filename",'s',false);
  imageASCIIOutputCmd->SetParameter(t10);

  G4UIparameter* t11 = new G4UIparameter("detector name",'s',false);
  imageASCIIOutputCmd->SetParameter(t11);

  mapCmd = new G4UIcommand("/scan/mapOn",this);
  mapCmd->SetGuidance("Generate a map of track ends .");
  mapCmd->SetGuidance("The map is generated by looping over the spherical coordinates theta and phi and");
  mapCmd->SetGuidance("recording the colour (visual attribute) of the logical volume where each track ends.");
  mapCmd->SetGuidance("Entering a negative value for nphi gives a polar plot of size ntheta x ntheta.");
  mapCmd->SetGuidance("WARNING: You may have to reset the angular distribution of your gps after this command!!!");
  mapCmd->AvailableForStates(G4State_Idle);
  G4UIparameter* m1 = new G4UIparameter("filename",'s',false);
  m1->SetDefaultValue("");
  m1->SetGuidance("The ouput file is a .ppm file");
  mapCmd->SetParameter(m1);
  G4UIparameter* m2 = new G4UIparameter("ntheta",'i',true);
  m2->SetDefaultValue("90");
  m2->SetGuidance("Number of steps in theta direction");
  mapCmd->SetParameter(m2);
  G4UIparameter* m3 = new G4UIparameter("nphi",'i',true);
  m3->SetDefaultValue("90");
  m3->SetGuidance("Number of steps in phi direction");
  mapCmd->SetParameter(m3);

  mapRotCmd = new G4UIcmdWith3Vector("/scan/mapAxisRotations",this);
  mapRotCmd->SetGuidance("Determine the orientation for the mapOn command");
  mapRotCmd->SetGuidance("Parameters: 3 angles (degrees): rotation around x,y and z axis (in that order)");
  mapRotCmd->SetParameterName("RX","RY","RZ",true,true);

  prefixCmd = new G4UIcmdWithAString("/scan/scanPrefix",this);
  prefixCmd->SetGuidance("Prefix for series of hit files produced for each point of the scan.");
  prefixCmd->SetGuidance("The output files will be named:   prefix###.hits");
  prefixCmd->SetParameterName("prefix",true);
  prefixCmd->SetDefaultValue("");


}

SEMScanMessenger::~SEMScanMessenger()
{
  delete scanCmd;
  delete imageCmd;
  delete scanXdirCmd;
  delete scanYdirCmd;
  delete scanNxCmd;
  delete scanNyCmd;
  delete scanLxCmd;
  delete scanLyCmd;
  delete scanOriginCmd;
  delete scanOutputCmd;
  delete imageOutputCmd;
  delete imageASCIIOutputCmd;
  delete mapCmd;
  delete mapRotCmd;
  delete prefixCmd;
}

void SEMScanMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==mapCmd ) {
     G4String filename;
     G4int    ntheta,nphi;
     const char* nv = (const char*)newValue;
     std::istringstream is((char*)nv);
     is >> filename >> ntheta >> nphi;
     scanManager->MapOn(filename,ntheta,nphi);
  }
  if( command==scanCmd )
  { 
     G4int nev;
     G4String filename;
     const char* nv = (const char*)newValue;
     std::istringstream is((char*)nv);
     is >> nev >> filename;
     if (filename == "") {
        G4cerr << "Performing normal linescan" << G4endl;
        scanManager->ScanOn(nev);
     } else {
        G4cerr << "Performing file-based linescan" << G4endl;
        scanManager->ScanOn(nev,filename);
     }
  }
  else if( command==imageCmd )
  {
     G4int nev;
     G4String filename;
     const char* nv = (const char*)newValue;
     std::istringstream is((char*)nv);
     is >> nev >> filename;
     G4cerr << "Intermediate output sent to file " << filename << G4endl;
     scanManager->ImageOn(nev,filename);
  }
  else if ( command==scanXdirCmd )
  { scanManager->SetScanXdir(scanXdirCmd->GetNew3VectorValue(newValue)); }
  else if ( command==mapRotCmd )
  { scanManager->SetMapRotation(mapRotCmd->GetNew3VectorValue(newValue)); }
  else if ( command==scanYdirCmd )
  { scanManager->SetScanYdir(scanYdirCmd->GetNew3VectorValue(newValue)); }
  else if ( command==scanNxCmd )
  { scanManager->SetScanNx(scanNxCmd->GetNewIntValue(newValue)); }
  else if ( command==scanNyCmd )
  { scanManager->SetScanNy(scanNyCmd->GetNewIntValue(newValue)); }
  else if ( command==scanLxCmd )
  { scanManager->SetScanLx(scanLxCmd->GetNewDoubleValue(newValue)); }
  else if ( command==scanLyCmd )
  { scanManager->SetScanLy(scanLyCmd->GetNewDoubleValue(newValue)); }
  else if ( command==scanOriginCmd )
  { scanManager->SetScanOrigin(scanOriginCmd->GetNew3VectorValue(newValue)); }
  else if ( command==prefixCmd )
  { scanManager->SetScanPrefix(newValue); }
  else if ( command==scanOutputCmd) { 
      string filename;
      G4String detectorname;
      const char* nv = (const char*) newValue;
      std::istringstream is(nv);
      is >> filename >> detectorname;
      scanManager->OutputLinescan(filename,detectorname);
  } else if ( command==imageOutputCmd) { 
      string filename;
      G4String detectorname;
      G4bool   autoscale;
      G4int    category;
      G4double min,max;
      const char* nv = (const char*) newValue;
      std::istringstream is(nv);
      is >> filename >> detectorname >> category >> autoscale >> min >> max;

      if (min>max) {
         G4double tmp;
         tmp = max;
         max = min;
         min = tmp;
         G4cout << "WARNING: min > max, values swapped" << G4endl;
      }
      
      scanManager->OutputPgm(filename,detectorname,category,autoscale,min,max);
  } else if ( command==imageASCIIOutputCmd) { 
      string filename;
      G4String detectorname;
      const char* nv = (const char*) newValue;
      std::istringstream is(nv);
      is >> filename >> detectorname;
      
      scanManager->OutputASCII(filename,detectorname);
  }
}

G4String SEMScanMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==scanXdirCmd )
  { cv = scanXdirCmd->ConvertToString(scanManager->GetScanXdir()); }
  else if( command==scanYdirCmd )
  { cv = scanYdirCmd->ConvertToString(scanManager->GetScanYdir()); }
  else if( command==scanNxCmd )
  { cv = scanNxCmd->ConvertToString(scanManager->GetScanNx()); }
  else if( command==scanNyCmd )
  { cv = scanNyCmd->ConvertToString(scanManager->GetScanNy()); }
  else if( command==scanLxCmd )
  { cv = scanLxCmd->ConvertToString(scanManager->GetScanLx()); }
  else if( command==scanLyCmd )
  { cv = scanLyCmd->ConvertToString(scanManager->GetScanLy()); }
  else if( command==scanOriginCmd )
  { cv = scanOriginCmd->ConvertToString(scanManager->GetScanOrigin()); }
  else if( command==prefixCmd )
  { cv = scanManager->GetScanPrefix(); }
  
  return cv;
}

