// EB
// EK
// 2007-01-04: Created special version that adds physical processes for newly created ions 'on the fly'.
// 2007-07-18: Added this history.
// See comments marked 'EK' further on in this file.
//
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       SEMPrimaryGeneratorActionMessenger.cc
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
#include "SEMPrimaryGeneratorAction.hh"
#include "SEMPrimaryGeneratorActionMessenger.hh"

#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4IonTable.hh"

#include "G4UIdirectory.hh"
#include "G4Tokenizer.hh"
//#include "G4SingleParticleSource.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ProcessManager.hh"
#include "G4Transportation.hh"
#include "CADPhysicsIonKill.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

///////////////////////////////////////////////////////////////////////////////
//
SEMPrimaryGeneratorActionMessenger::SEMPrimaryGeneratorActionMessenger(SEMPrimaryGeneratorAction *fPGA) 
    : target(fPGA),fShootIon(false)
{
  fParticleGun = target->GetGpsGun()->GetCurrentSource(); // Returns the gps gun (needed for ion command)
  //particleTable = G4ParticleTable::GetParticleTable();
  ionTable = G4IonTable::GetIonTable();

  // Overwrite original gps/ion command with our new version
  ionCmd = new G4UIcommand("/gps/ion2",this);
  ionCmd->SetGuidance("Set properties of ion to be generated (version including Transportation process).");
  ionCmd->SetGuidance("[usage] /gps/ion Z A Q E");
  ionCmd->SetGuidance("        Z:(int) AtomicNumber");
  ionCmd->SetGuidance("        A:(int) AtomicMass");
  ionCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
  ionCmd->SetGuidance("        E:(double) Excitation energy (in keV)");
  
  G4UIparameter* param;
  param = new G4UIparameter("Z",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("A",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("Q",'i',true);
  param->SetDefaultValue("0");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  ionCmd->SetParameter(param);

  gpsgunCmd = new G4UIcmdWithABool("/gps/gpsgun",this);
  gpsgunCmd->SetGuidance("Select gps gun (default) or special gun which reads input from a file generated with OutputHits");
  gpsgunCmd->SetParameterName("gpsgun",true);
  gpsgunCmd->SetDefaultValue(true);

  hitsfileCmd = new G4UIcmdWithAString("/gps/hitsfile",this);
  hitsfileCmd->SetGuidance("Determines which file the primary particles are read from (if gpsgun is set to false)");
  hitsfileCmd->SetGuidance("Enter a filename as parameter");
  hitsfileCmd->SetParameterName("filename",false);

}

SEMPrimaryGeneratorActionMessenger::~SEMPrimaryGeneratorActionMessenger()
{
  delete ionCmd;
  delete gpsgunCmd;
  delete hitsfileCmd;
}

void SEMPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
  if (command==ionCmd) { IonCommand(newValues); }
  if (command==gpsgunCmd) { target->SetGun(gpsgunCmd->GetNewBoolValue(newValues)); }
  if (command==hitsfileCmd) { target->ReadHitsFile(newValues); }
}

G4String SEMPrimaryGeneratorActionMessenger::GetCurrentValue(G4UIcommand *)
{
  G4String cv;
  
  cv = "Not implemented yet";

  return cv;
}

void SEMPrimaryGeneratorActionMessenger::IonCommand(G4String newValues)
{
  fShootIon = true;

  if (fShootIon)
  {
    G4Tokenizer next( newValues );
    // check argument
    fAtomicNumber = StoI(next());
    fAtomicMass = StoI(next());
    G4String sQ = next();
    if (sQ.isNull())
    {
       fIonCharge = fAtomicNumber;
    }
    else
    {
       fIonCharge = StoI(sQ);
       G4cout << fIonCharge << G4endl;
       sQ = next();
       if (sQ.isNull())
          {
             fIonExciteEnergy = 0.0;
          }
          else
          {
             fIonExciteEnergy = StoD(sQ) * keV;
          }
    }
    G4ParticleDefinition* ion;
    ion =  ionTable->GetIon( fAtomicNumber, fAtomicMass, fIonExciteEnergy);
    if (ion==0)
    {
      G4cout << "Ion with Z=" << fAtomicNumber;
      G4cout << " A=" << fAtomicMass << "is not be defined" << G4endl;    
    }
    else
    {
      G4cerr << "Using extended Geant4 Ion definition (including Transportation)" << G4endl;
      // EK, this is the new part. When a particular ion is selected as the primary particle,
      // it is checked whether this ion already has processes assigned to it. If not, 
      // new Transportation and CADPhysicsIonKill processes are created and assigned to the particle
      // 'on the fly'.
      G4ProcessManager* pmanager = ion->GetProcessManager();
      if ( pmanager && pmanager->GetProcessListLength()==0) {
         G4Transportation* theTransportationProcess = new G4Transportation();
         G4VProcess* ionkill                        = new CADPhysicsIonKill();
         pmanager ->AddProcess(theTransportationProcess);
         pmanager ->AddDiscreteProcess(ionkill);
         pmanager ->SetProcessOrdering(theTransportationProcess, idxAlongStep,0);
         pmanager ->SetProcessOrdering(theTransportationProcess, idxPostStep,0);
         pmanager ->SetProcessOrdering(ionkill,                  idxPostStep,1);
      }
      // EK, End of the new part
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(fIonCharge*eplus);
    }
  }
  else
  {
    G4cout << "Set /gps/particle to ion before using /gps/ion command";
    G4cout << G4endl; 
  }
}
