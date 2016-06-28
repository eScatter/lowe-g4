///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       SEMPrimaryGeneratorActionMessenger.hh
//
//
///////////////////////////////////////////////////////////////////////////////
//
//
#ifndef SEMPrimaryGeneratorActionMessenger_h
#define SEMPrimaryGeneratorActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"
#include "CADPhysicsUnits.hh"

class G4ParticleTable;
class G4IonTable;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class SEMPrimaryGeneratorAction;
//class G4GeneralParticleSource;
class G4SingleParticleSource;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class SEMPrimaryGeneratorActionMessenger: public G4UImessenger
{
public:
  SEMPrimaryGeneratorActionMessenger(SEMPrimaryGeneratorAction*);
  ~SEMPrimaryGeneratorActionMessenger();

  void SetNewValue(G4UIcommand *command, G4String newValues);
  //    Identifies the command which has been invoked by the user, extracts the
  //    parameters associated with that command (held in newValues), and uses
  //    these values with the appropriate member function of G4GeneralParticleSource.
  G4String GetCurrentValue(G4UIcommand *command);

private:
  void IonCommand(G4String newValues);

private:
//  G4GeneralParticleSource *fGPS;
  SEMPrimaryGeneratorAction *target;
  G4SingleParticleSource  *fParticleGun;
  G4ParticleTable *particleTable;
  G4IonTable *ionTable;

    
private: //commands
  G4UIcommand         *ionCmd;
  G4UIcmdWithABool    *gpsgunCmd;
  G4UIcmdWithAString  *hitsfileCmd;

private: // for ion shooting
  G4bool   fShootIon; 
  G4int    fAtomicNumber;
  G4int    fAtomicMass;
  G4int    fIonCharge;
  G4double fIonExciteEnergy;

};

#endif

