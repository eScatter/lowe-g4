#ifndef SEMPrimaryGeneratorAction_h
#define SEMPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include <vector>

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

using namespace std;

class G4Event;
class G4SingleParticleSource;
class SEMPrimaryGeneratorActionMessenger;

class SEMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    SEMPrimaryGeneratorAction();    
   ~SEMPrimaryGeneratorAction();

  public:
    G4ThreeVector position;
    void GeneratePrimaries(G4Event*);
    void SetPosition(G4ThreeVector);
    G4ThreeVector GetPosition();
    G4double GetParticleEnergy();
    G4int GetNumberOfParticles();
    void SetNumberOfParticles(G4int n);
    G4ThreeVector GetParticleMomentumDirection();
    G4SingleParticleSource* GetCurrentSource();

    inline void SetGun(G4bool g) { G4cerr << "Value of gpsgun is now " << gpsgun << G4endl; gpsgun=g; G4cerr << "Setting gpsgun to " << gpsgun << G4endl;}
    inline G4bool GetGun() { return gpsgun;}
    void ReadHitsFile(G4String filename);

    inline G4GeneralParticleSource* GetGpsGun() { return particleGun;}
    inline G4ParticleGun* GetHitsGun() { return particleGun2;}

protected:
    G4GeneralParticleSource* particleGun;
    G4ParticleGun* particleGun2;
    SEMPrimaryGeneratorActionMessenger* messenger;
    G4bool gpsgun; // Determines which gun to use
    G4String hitsfilename;

    vector <G4double> KinEnergy;
    vector <G4ThreeVector> Position;
    vector <G4ThreeVector> Direction;
    vector <G4String> ParticleType;

};

#endif


