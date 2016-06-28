#include "SEMDetectorConstruction.hh"
#include "CADPhysicsPhysicsList.hh"
#include "SEMPrimaryGeneratorAction.hh"
#include "SEMRunAction.hh"
#include "SEMTrackingAction.hh"
#include "SEMEventAction.hh"
#include "SEMSteppingAction.hh"
#include "SEMSteppingVerbose.hh"
#include "G4RunManager.hh"

//#include "SEMScanManager.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
// This visualisation executive takes care of transparency in the model. Also it allows visualization on nm scale
#include "SEMVisExecutive.hh"
#endif

int main(int argc,char** argv) {

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SEMSteppingVerbose);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  SEMDetectorConstruction* SEMdetector;
  if (argc > 1) {
     G4String fileName=argv[1];
     G4cout << "GDML file name = " << fileName << G4endl;
     SEMdetector = new SEMDetectorConstruction(fileName);
  } else {
     SEMdetector = new SEMDetectorConstruction;
  }
  runManager->SetUserInitialization(SEMdetector);
  runManager->SetUserInitialization(new CADPhysicsPhysicsList);
  
  //Initialize G4 kernel
  runManager->Initialize();
      
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new SEMVisExecutive;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new SEMPrimaryGeneratorAction);
  runManager->SetUserAction(new SEMRunAction);  
  runManager->SetUserAction(new SEMTrackingAction);
  runManager->SetUserAction(new SEMEventAction);
  runManager->SetUserAction(new SEMSteppingAction);
//  SEMScanManager* scanmanager = new SEMScanManager;

  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  // Define (G)UI terminal for interactive mode  
  // Usage:
  // SEM [gdmlfile [macfile]]
  // Behaviour:
  // if no arguments are given (argc==1) : execute SEM.mac using SEM.gdml and end with a user prompt
  // if one argument is given (argc==2): execute SEM.mac using given gdml file and end with a user prompt
  // if two arguments are given (argc==3): execute given macfile uwing given gdml file and exit 
  if(argc<=2) { 
     // G4UIterminal is a (dumb) terminal.
     G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
     G4cout << "Using G4UItcsh" << G4endl;
     session = new G4UIterminal(new G4UItcsh);      
#else
     G4cout << "Using G4UIterminal" << G4endl;
     session = new G4UIterminal();
#endif    

     UI->ApplyCommand("/control/execute SEM.mac");    
     session->SessionStart();
     delete session;
  } else { 
     // Batch mode
     G4String command = "/control/execute ";
     G4String fileName = argv[2];
     UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

